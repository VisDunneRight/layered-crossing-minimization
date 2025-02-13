import csv


def analyze_results():
    with open("results/_results.csv") as fd:
        rdr = csv.reader(fd)
        next(rdr)
        sum1, sum2 = 0, 0
        diff = 0
        diff_div = 0
        ct = 0
        for ln in rdr:
            if float(ln[2]) > 300 or float(ln[3]) > 300:
                print(ln)
                continue
            sum1 += float(ln[2])
            sum2 += float(ln[3])
            diff += float(ln[2]) - float(ln[3])
            diff_div += (float(ln[2]) - float(ln[3])) / float(ln[2])
            ct += 1
        print(sum1, sum2, ct)
        print("Avg1:", round(sum1 / ct, 5))
        print("Avg2:", round(sum2 / ct, 5))
        print("Avg T1-T2:", round(diff / ct, 5))
        print("Avg (T1-T2)/T1:", round(diff_div / ct, 5))
        print(f"\tas percent: {round(diff_div / ct, 5) * 100}%")


if __name__ == '__main__':
    analyze_results()
