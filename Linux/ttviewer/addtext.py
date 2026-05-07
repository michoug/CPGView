import sys
from pymupdf import pymupdf
from Bio import SeqIO

arch1 = sys.argv[1]
name = sys.argv[1]

for rec in SeqIO.parse(arch1,'gb'):
    m = (rec.annotations)["organism"]
    
name1 = str(str(name).replace(".gbf",""))

doc = pymupdf.open(name1 + '_circle0.pdf') 

def add_footer(pdf):
    for i in range(0, pdf.page_count):
        page = pdf[i]
        page.insert_text((15,50),m,fontname="Times-Italic",fontsize=15)
        #page.insertText((15,50),m,fontname="Times-Italic",fontsize=15)


add_footer(doc)
result = name1 + '_circle.pdf'
doc.save(result)
