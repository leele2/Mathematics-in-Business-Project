# import module
from pdf2image import convert_from_path
from pathlib import Path
from os.path import exists
from os import mkdir

# Directories
cwd = str(Path(__file__).parent.absolute())
pdf_file = cwd + "\\main.pdf"
output_dir = cwd + "\\Images\\"
readme_dir = "\\".join(cwd.split("\\")[:-3]) # -3 represents up 4 directories
 
 
# Store Pdf with convert_from_path function
images = convert_from_path(pdf_file)
txt_out = []
txt_out.append("#" + cwd.split("\\")[-4] + "\n")
# txt_out.append("Statistical Modelling in Space and Time CW2\n")

# Create folder to store images if not already created
if not exists(output_dir):
    mkdir(output_dir)

for i in range(len(images)):
    #Save pages as images in the pdf
    images[i].save(output_dir + 'page'+ str(i) +'.png', 'PNG')
    #String for README.md
    txt_out.append("![page" + str(i) + "](Latex_Files/Main/z_output/Images/page" + str(i) + ".png)")
    txt_out.append("***")

with open(readme_dir + "/README.md", "w") as output:
    output.write("\n".join(txt_out))