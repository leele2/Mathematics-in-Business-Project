# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 18:05:55 2022
Scraping and plotting data from theocc.com
@author: leele2
"""
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from bs4 import BeautifulSoup
from os import listdir
from pathlib import Path

def string_to_int(string):
    if not string:
        return 0
    if string[0] == "$":
        string = string[1:]
    return int(string.replace(',', ''))

# Ensuring html file is on path
"""
html file taken from:
https://www.theocc.com/Market-Data/Market-Data-Reports/Volume-and-Open-Interest/Historical-Volume-Statistics
"""
html_file = [f for f in listdir('.') if f.endswith('.html')]
if len(html_file) != 1:
    raise ValueError('should be exactly one html file in the current directory')
# Importing html file
with open(html_file[0], 'r', encoding='utf-8') as file:
    soup = BeautifulSoup(file, 'lxml')
# Finding Table
table = soup.find_all("table")
# Selecting Data
table_data = table[0].find_all("tr")
# Finding Table Names
table_names = []
for names in table_data[0]:
    table_names.append(str(names.find("span").string))
del table_names[0]
del table_names[-1]
# Finding Headers
headers = []
for header in table_data[1]:
    headers.append(str(header.find("span").string))
# Populating Data
Data = {}
for rows in table_data[2:]:
    temp = {}
    for i, data in enumerate(rows):
        if i > 6 and i < 10 or i == 0:
            temp[headers[i]] = data.get_text()
    Date = datetime.strptime(temp["Date"], "%Y")
    del temp["Date"]
    Data[Date] = {k: string_to_int(v) for k, v in temp.items()}
# Converting data into pandas dataframes
Data = pd.DataFrame.from_dict(Data, orient="index")
# Plotting time-series
plt.rc('font', family='serif')
fig, ax1 = plt.subplots(figsize=(6, 4.5))
ax1.plot(Data['Options']/1e6, linewidth=2, zorder=1)
ax1.set_title('Average Option Volume (millions)')
ax1.set_xlabel('Date')
ax1.xaxis.set_minor_formatter(mdates.DateFormatter("%Y"))
ax1.set_ylabel('Average Daily Volume')
ax1.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}M'))
fig.autofmt_xdate()
plt.tight_layout()
plt.show()
# Saving plot
cwd = str(Path(__file__).parent.parent.absolute())
fig.savefig(cwd + "\Latex_Files\Main\Chapters\C1\plots\OptionVolume.png")
