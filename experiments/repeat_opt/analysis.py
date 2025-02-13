import csv


def analyze_results():
    with open("results/_results.csv") as fd:
        rdr = csv.reader(fd)
        next(rdr)
        sum1, sum2, sum3 = [0] * 8, [0] * 8, [0] * 8
        fms, smt, tmf = [0] * 8, [0] * 8, [0] * 8
        fmsd, smtd, tmfd = [0] * 8, [0] * 8, [0] * 8
        ct = [0] * 8
        idx_map = {}
        next_int = 0
        for ln in rdr:
            if ln[1][-2:] not in idx_map:
                idx_map[ln[1][-2:]] = next_int
                next_int += 1
            sum1[idx_map[ln[1][-2:]]] += float(ln[2])
            sum2[idx_map[ln[1][-2:]]] += float(ln[3])
            sum3[idx_map[ln[1][-2:]]] += float(ln[4])
            fms[idx_map[ln[1][-2:]]] += float(ln[2]) - float(ln[3])
            smt[idx_map[ln[1][-2:]]] += float(ln[3]) - float(ln[4])
            tmf[idx_map[ln[1][-2:]]] += float(ln[4]) - float(ln[2])
            fmsd[idx_map[ln[1][-2:]]] += (float(ln[2]) - float(ln[3])) / float(ln[2])
            smtd[idx_map[ln[1][-2:]]] += (float(ln[3]) - float(ln[4])) / float(ln[3])
            tmfd[idx_map[ln[1][-2:]]] += (float(ln[4]) - float(ln[2])) / float(ln[4])
            ct[idx_map[ln[1][-2:]]] += 1
        print(sum1, sum2, sum3, ct)
        print("Avg1:", [round(sum1[i] / ct[i], 3) for i in range(8)])
        print("Avg2:", [round(sum2[i] / ct[i], 3) for i in range(8)])
        print("Avg3:", [round(sum3[i] / ct[i], 3) for i in range(8)])
        print("Avg T1-T2:", [round(fms[i] / ct[i], 3) for i in range(8)])
        print("Avg T2-T3:", [round(smt[i] / ct[i], 3) for i in range(8)])
        print("Avg T3-T1:", [round(tmf[i] / ct[i], 3) for i in range(8)])
        print("Avg (T1-T2)/T1:", [round(fmsd[i] / ct[i], 3) for i in range(8)])
        print("Avg (T2-T3)/T2:", [round(smtd[i] / ct[i], 3) for i in range(8)])
        print("Avg (T3-T1)/T3:", [round(tmfd[i] / ct[i], 3) for i in range(8)])


def delete_bad_rows():
    with open("results/_results.csv") as fd:
        lines = list(fd.readlines())
        for i in range(len(lines)-1, 1, -1):
            spline = lines[i].removesuffix('\n').split(',')
            print(spline)
            if float(spline[2]) >= 60 or float(spline[3]) >= 60 or float(spline[4]) >= 60:
                lines.pop(i)
        with open("results/_results.csv", 'w') as fd2:
            fd2.writelines(lines)


def reindex():
    with open("results/_results.csv") as fd:
        lines = list(fd.readlines())
        last_idx = 0
        for i in range(1, len(lines)):
            spline = lines[i].split(',')
            idx = int(spline[0])
            if idx != last_idx + 1:
                spline[0] = str(last_idx + 1)
                lines[i] = ','.join(spline)
            last_idx += 1
        with open("results/_results.csv", 'w') as fd2:
            fd2.writelines(lines)


if __name__ == '__main__':
    analyze_results()
