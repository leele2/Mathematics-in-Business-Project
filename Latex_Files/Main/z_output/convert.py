# import module
from pdf2image import convert_from_path
import sys

pdf_file = sys.argv[1] 
output_dir = sys.argv[2] 
 
 
# Store Pdf with convert_from_path function
images = convert_from_path(pdf_file)
 
for i in range(len(images)):
   
      # Save pages as images in the pdf
    images[i].save(output_dir + 'page'+ str(i) +'.jpg', 'JPEG')