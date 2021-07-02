import openpyxl
wb = openpyxl.load_workbook('data/src/sample.xlsx')
sheet = wb['sheet1']
kai = []
for i in range (2,1000):
    kai.append(sheet.value(row=2,column=2))
