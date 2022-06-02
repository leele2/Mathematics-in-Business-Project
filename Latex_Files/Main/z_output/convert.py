# import module
from pdf2image import convert_from_path
import sys

pdf_file = sys.argv[1]
output_dir = sys.argv[2]
readme_dir = sys.argv[3]
 
 
# Store Pdf with convert_from_path function
images = convert_from_path(pdf_file)
txt_out = []
txt_out.append("# Mathematics-in-Business-Project\n")
 
for i in range(len(images)):
    #Save pages as images in the pdf
    images[i].save(output_dir + 'page'+ str(i) +'.jpg', 'JPEG')
    #String for README.md
    txt_out.append("![img](Latex_Files/Main/z_output/Images/page" + i + ".jpg)")

