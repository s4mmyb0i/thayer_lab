start_res=131 #beginning index for 1nf3 
end_res = 253 # ending index for 1nf3
residuemap = {} 

# binding pocket only map: 
p_pocket_list = [3, 4, 6, 9, 10, 11, 13, 14, 17, 19, 15, 18, 26, 78]
pocketmap = {} 
for p_num in p_pocket_list:
    ref_num = p_num + 130
    pocketmap[p_num] = ref_num

print(f"{'P.pdb Pocket #'} -> {'1NF3 Pocket #'}")
print("-" * 34)
for p_num, ref_num in pocketmap.items():
    print(f"{p_num} -> {ref_num}")

# entire map: 
for ref_res in range (start_res, end_res +1 ): #including 253 as part of range 
    p_res = ref_res - 130 ##offset by 130 
    if p_res > 0: 
        residuemap [p_res] = ref_res

""" uncomment if u want to see the entire dictionary b/t p and 1nf3 indices 
print(f"{'P.pdb res'} -> {'1NF3 Number'}")
print("-" * 28)
for p_res, ref_res in sorted(residuemap.items()):
    print(f"{p_res} -> {ref_res}")

""" 

while True: 
    user_input = input("enter P indices separated by comma or space or 'q' to quit: ")
    if user_input.lower() == 'q':
        break
    try: 
        cleanedInput = user_input.replace (',', ' ')
        inputs = [int(x) for x in cleanedInput.split()]
    except ValueError: 
        print ("Pls check if all inputs are integers \n")
        continue 
    outputs = {} 
    missingRes = [] 
    for pRes_key in inputs: 
        if pRes_key in residuemap: 
            outputs[pRes_key] = residuemap[pRes_key] 
        else: 
            missingRes.append(pRes_key)
    print (outputs)
    if missingRes: 
        print (f"residues {missingRes} are not part of the map.")
    print ("-" * 5 + "\n")

