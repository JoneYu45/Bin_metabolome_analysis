#Environment setup
import numpy as np
import re
import pandas as pd
from optparse import OptionParser

def main(input, output):
    # Bin info input
    bin_ko_info = np.empty(shape=(1, 7))
    with open(input) as input_info:
        print('Inputting bin data...')
        for bin_info in input_info:
            new_info = bin_info.split(sep='\t')
            if len(new_info) > 1:
                new_info = [new_info[0].strip(), new_info[1].strip(), '__', '__', '__', '__', '__']
            else:
                new_info = [new_info[0].strip(), '__', '__', '__', '__', '__', '__']
            new_info = np.array(new_info)
            bin_ko_info = np.row_stack((bin_ko_info, new_info))
        print('Done.')
    bin_ko_info = bin_ko_info[1:, :]

    # Save temp
    bin_ko_info = pd.DataFrame(bin_ko_info)
    bin_ko_info.to_csv('../outputs/temp.csv')
    bin_ko_info = pd.read_csv('../outputs/temp.csv', index_col=0)
    bin_ko_info = np.array(bin_ko_info)

    # Map EC number according to the KO number against KEGG
    print('Mapping EC number according to the KO number against KEGG...')
    with open('../references/ko00001.keg') as ko:
        n = 0
        for ko_info in ko:
            # Locate EC number
            if re.search('D      K[0-9]+', ko_info) != None:
                for i in np.where(bin_ko_info[:, 1] != '__')[0]:
                    if bin_ko_info[i, 1] in ko_info:
                        EC_loc = re.search('\[EC:.+\]', ko_info)
                        if EC_loc != None:
                            bin_ko_info[i, 2] = EC_loc[0]

            # Report progress
            n = n + 1
            print('\rMapping %.2f %%' % float(n / 54096 * 100), end='', flush=True)

    # Save temp
    bin_ko_info = pd.DataFrame(bin_ko_info)
    bin_ko_info.to_csv('../outputs/temp.csv')
    bin_ko_info = pd.read_csv('../outputs/temp.csv', index_col=0)
    bin_ko_info = np.array(bin_ko_info)

    # Map compound number according to EC number
    print('\nMapping compound number according to EC number...')
    with open('../references/br08202.keg') as br:
        n = 1
        target_EC = 0
        for br_info in br:
            # Find target EC number
            if re.search('B  [0-9].', br_info) != None:
                target_EC = 0
                for i in np.where(bin_ko_info[:, 2] != '__')[0]:
                    EC_number = '\.'.join(re.split('\.', br_info.strip().replace('B  ', '')))
                    if re.search(EC_number + '[^0-9]', bin_ko_info[i, 2]) != None:
                        target_EC = 1
                        break

            # Map the compound ID under the target EC
            if target_EC == 1:
                C_numbers = re.findall('C[0-9][0-9][0-9][0-9][0-9]', br_info)
                if len(C_numbers) > 2:
                    bin_ko_info[i, 3] = '_'.join((bin_ko_info[i, 3], C_numbers[0], C_numbers[2]))
                else:
                    if len(C_numbers) > 0:
                        bin_ko_info[i, 3] = '_'.join((bin_ko_info[i, 3], C_numbers[0]))

            # Report progress
            n = n + 1
            print('\rMapping %.2f %%' % float(n / 1790 * 100), end='', flush=True)

    #Save temp
    bin_ko_info = pd.DataFrame(bin_ko_info)
    bin_ko_info.to_csv('../outputs/temp.csv')
    bin_ko_info = pd.read_csv('../outputs/temp.csv', index_col=0)
    bin_ko_info = np.array(bin_ko_info)

    # Map compound numbers against HMDB to gain formula and chemical taxonomy
    print('\nMapping compound numbers against HMDB to gain formula and chemical taxonomy...')
    with open('../references/hmdb_metabolites.xml', encoding='utf-8') as hmdb:
        print('Preparing database...')
        length = 0
        for hmdb_info in hmdb:
            length = length + 1

    with open('../references/hmdb_metabolites.xml', encoding='utf-8') as hmdb:
        print('Start mapping...')
        temp_info = []
        n = 0
        for hmdb_info in hmdb:
            # Collect wanted info of a metabolite
            if re.search('  <chemical_formula>', hmdb_info) != None:
                temp_info.append(hmdb_info.strip().replace('<chemical_formula>', '').replace('</chemical_formula>', ''))
            if re.search('    <direct_parent>', hmdb_info) != None:
                temp_info.append(hmdb_info.strip().replace('<direct_parent>', '').replace('</direct_parent>', ''))
            if re.search('    <super_class>', hmdb_info) != None:
                temp_info.append(hmdb_info.strip().replace('<super_class>', '').replace('</super_class>', ''))

            # See such info useful or not
            if re.search('  <kegg_id>', hmdb_info) != None:
                current_C_num = hmdb_info.strip().replace('<kegg_id>', '').replace('</kegg_id>', '')
                for i in np.where(bin_ko_info[:, 3] != '__')[0]:
                    all_C_num = bin_ko_info[i, 3].split(sep='_')[3:]
                    if current_C_num in all_C_num:
                        if len(temp_info) < 3:
                            pass
                        else:
                            bin_ko_info[i, 4] = '_'.join((bin_ko_info[i, 4], temp_info[0]))
                            bin_ko_info[i, 5] = '_'.join((bin_ko_info[i, 5], temp_info[1]))
                            bin_ko_info[i, 6] = '_'.join((bin_ko_info[i, 6], temp_info[2]))
            if re.search('</metabolite>', hmdb_info) != None:
                temp_info = []

            # Report progress
            n = n + 1
            print('\rMapping %.2f %%' % float(n / length * 100), end='', flush=True)
    print('\nDone.', flush=True)

    # Temp_file output
    bin_ko_info = pd.DataFrame(bin_ko_info)
    bin_ko_info.to_csv(output)

if __name__ == '__main__':
    # Imnput data and setting
    parse = OptionParser()
    parse.add_option('-I', '--input', dest='input', default='../inputs/6462_u_ko.txt')
    parse.add_option('-O', '--output', dest='output', default='../outputs/temp.csv')
    (options, args) = parse.parse_args()

    input = options.input
    output = options.output

    main(input, output)