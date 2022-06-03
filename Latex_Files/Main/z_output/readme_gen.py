# import module
from pdf2image import convert_from_path
import sys

pdf_file = sys.argv[1]
output_dir = sys.argv[2]
readme_dir = "C:\\Users\\dj-lu\OneDrive - University of Exeter\\University of Exeter\\05 - Fifth Year\\Mathematics in Business Project\\"
 
 
# Store Pdf with convert_from_path function
images = convert_from_path(pdf_file)
txt_out = []
txt_out.append("# Mathematics-in-Business-Project\n")
 
for i in range(len(images)):
    #Save pages as images in the pdf
    images[i].save(output_dir + 'page'+ str(i) +'.png', 'PNG')
    #String for README.md
    txt_out.append("![page" + str(i) + "](Latex_Files/Main/z_output/Images/page" + str(i) + ".png)")
    txt_out.append("***")

with open(readme_dir + "README.md", "w") as output:
    output.write("\n".join(txt_out))