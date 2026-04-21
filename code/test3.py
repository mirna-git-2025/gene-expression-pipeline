cell_count = int(input("Enter cell count: "))
print(f"You entered {cell_count} cells")

def doublecell(countc):
    return countc * 2 

if  cell_count > 100 :
    print("big problem")
else:
    print("ok no problem")


print(" formule " , doublecell(cell_count))    


