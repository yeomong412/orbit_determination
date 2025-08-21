input_filename = 'data/nasa_raw_input.txt'
output_lines = []

with open(input_filename, 'r') as f:
    lines = f.readlines()

# '$$SOE'와 '$$EOE' 사이의 데이터만 사용
data_lines = []
soe_found = False
for line in lines:
    if '$$SOE' in line:
        soe_found = True
        continue
    if '$$EOE' in line:
        break
    if soe_found:
        data_lines.append(line.strip())

count = 0

for line in data_lines:
    if not line:
        continue
    # 공백을 기준으로 분리하면 여러 개의 공백이 있는 경우에도 잘 분리됨
    parts = line.split()
    jd = parts[0]
    ra = parts[1]
    dec = parts[2]
    delta = parts[3]
    # 원하는 순서: R.A., DEC, delta, JDUT
    output_lines.append(f"{ra} {dec} {delta} {jd}")
    count += 1

inputf = open('data\input.txt', 'w')

for line in output_lines:
    print(line)
    inputf.write(line)
    inputf.write('\n')

inputf.close()

print(f"총 {count}개")