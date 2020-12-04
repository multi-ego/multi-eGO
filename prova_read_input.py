
with open('GRETA/native/topol.top') as input_data:
    # Skips text before the beginning of the interesting block:
    for line in input_data:
        if line.strip() == '[ atoms ]':
            break
    # Reads text until the end of the block:
    for line in input_data:
        if line.strip() == '':
            break
        print(line)