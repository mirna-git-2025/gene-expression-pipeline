cell_count = 5000
temperature = 37.5
gene_name = "BRCA1"
print(f"Cell count: {cell_count}, Temp: {temperature}, Gene: {gene_name}")
culture_volume = 2.5  # liters
nutrient_concentration = 0.1  # g/L
print(f"Culture: {culture_volume} L, Nutrient: {nutrient_concentration} g/L")

population = 100
print("Population growth:")
while population < 10000: 
    population *= 2
    print(population)

bacteria = 50
hours = 0
while bacteria < 1000:
    bacteria *= 2
    hours += 1
    print(f"Hour {hours}: {bacteria} bacteria")


name = "helene"
age=20
print(" researcher " , name , "age" , age)
print(f" researcher  {name} : her age {age}")
genes = ["G1", "G2", "G3"]
for g in genes:
    print("Gene:", g)
# List of first 20 genes
genes = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10',
         'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18', 'G19', 'G20']

# Loop from i = 3 to i = 7
for i in range(3, 8):  # range is inclusive of 3, exclusive of 8
    if i < len(genes):  # make sure we don't go out of bounds
        print(f"Index {i}: {genes[i]}")
    
    

def calculate_concentration(mass, volume):
    return mass / volume

print("Concentration1:", calculate_concentration(10, 2))
print("Concentration2:", calculate_concentration(30, 5))
def od_per_cell(od, cells):
    return od / cells

print(od_per_cell(1.2, 1e6))

expression = 10
if expression >= 10:
    if expression > 20:
        print("Extremely high")
    else:
        print("Moderately high")

else:
    
    print("asmas high")

# List
enzymes = ["DNA polymerase", "RNA polymerase", "Ligase"]
print(enzymes[0])
print(enzymes[1])

# Tuple (immutable)
temperature_range = (20, 37, 42)
print(temperature_range)

# Dictionary (key-value)
gene_expression = {"G1": 50, "G2": 120, "G3": 75}
print(gene_expression["G2"])
print(gene_expression["G3"])
        
gene_expression = {"G1": 50, "G2": 120, "G3": 75}

for gene, value in gene_expression.items():
    print("gene" ,gene, "value " , value)
gene_expression = {"G1": 50, "G2": 120, "G3": 75}

for i, (gene, value) in enumerate(gene_expression.items()):
    if i > 2 :
        print(i, gene, value)


cell_count = int(input("Enter cell count: "))
print(f"You entered {cell_count} cells")
        
