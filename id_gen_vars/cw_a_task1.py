import sys
import csv

from scipy.stats.distributions import chi2

# Returns dict of SNPs and obseved alleles from tab seperated genotype file
def get_snp_dict(filename):
    with open(filename) as geno_file:
        geno_reader = csv.reader(geno_file, delimiter='	')
        # Get the names of the SNPs from the file header
        snps = next(geno_reader)
        # Get number of SNPs
        n_snps = len(snps)

        # Initialise SNP dict
        snp_dict = {}
        for snp in snps:
            snp_dict[snp] = []
        
        # Iterate through rows and populate SNP dict
        for row in geno_reader:
            for i in range(n_snps):
                snp_dict[snps[i]].append(row[i])
        
        return snp_dict

# Helper function which adds unoberved alleles necessary to calculate if in HWE
def get_obs_missing_alleles(obs_dict):
    bases = set()
    for allele in obs_dict.keys():
        for base in allele:
            bases.add(base)

    if len(bases) != 2:
        print(bases)
        raise Exception('Error in observed alleles, skipping')

    bases = list(bases)
    alleles = [f'{bases[0]}{bases[0]}',f'{bases[0]}{bases[1]}',f'{bases[1]}{bases[1]}']

    for allele in alleles:
        if allele not in obs_dict.keys():
            allele = allele[::-1]
            if allele not in obs_dict.keys():
                obs_dict[allele] = 0

    return obs_dict

# Returns a list of SNPs not in HWE from a genotype file
def get_nonhwe_snps(filename):
    # Initialise list which will store SNPs not in HWE
    nonhwe_snps = []
    # Get dict of SNPs and observed alleles
    snp_dict = get_snp_dict(filename)
        
    for snp in snp_dict:
        try:
            # Create a dict of observed allele numbers
            obs_dict = {}
            for ob in snp_dict[snp]:
                try:
                    if ob not in ['NA','na',None]: # Ignore NA observations
                        ob = ob.upper()
                        if ob not in obs_dict.keys():
                            ob = ob[::-1]
                            if ob not in obs_dict.keys():
                                obs_dict[ob] = 1
                            else:
                                obs_dict[ob] += 1
                except:
                    pass
            
            # Add the missing combinations as observerd at 0
            if len(obs_dict.keys()) < 3:
                obs_dict = get_obs_missing_alleles(obs_dict)

            # If there are more than 3 kinds of obsevation, ignore
            if len(obs_dict.keys()) > 3:
                raise Exception('Incorrect number of alleles')
            
            # Tally up the obervations
            obs_homoz_1 = None
            obs_heteroz = None
            obs_homoz_2 = None
            
            for allele in obs_dict.keys():
                if allele[0] != allele[1]:
                    obs_heteroz = obs_dict[allele]
                elif obs_homoz_1 == None:
                    obs_homoz_1 = obs_dict[allele]
                else:
                    obs_homoz_2 = obs_dict[allele]
           
            if None in (obs_homoz_1, obs_homoz_2, obs_heteroz):
                raise Exception('error calculating if in HWE')
            
            # Carry out a Chi Squared test to see if in HWE
            n_obs = obs_homoz_1 + obs_heteroz + obs_homoz_2

            # Calculate minor allele frequency
            ma_freq =  (2*obs_homoz_1 + obs_heteroz)/(2*n_obs)

            # Calculate expected observations from minor allele frequency
            exp_counts = [n_obs*(ma_freq**2), 2*n_obs*ma_freq*(1-ma_freq), n_obs*((1-ma_freq)**2)]

            obs_counts = [obs_homoz_1, obs_heteroz, obs_homoz_2]

            # Compute Chi squared statistic
            chi_sqrd_stat = 0
            for i in range(len(exp_counts)):
                obs_minus_exp = obs_counts[i] - exp_counts[i]
                chi_sqrd_stat += (obs_minus_exp**2)/exp_counts[i]
            
            # Calculate p-value from Chi Squared statistic and degrees of freedom (in this case 1)
            p_value = chi2.sf(chi_sqrd_stat,1)

            # Append snp to list of SNPs not in HWE
            if p_value >= 0.05:
                nonhwe_snps.append(snp)

        except Exception as err:
            # Ignore any SNPS that throw errors
            print(err)
            pass
    
    return nonhwe_snps

if __name__ == '__main__':
    filename = sys.argv[1]

    nonhwe_snps = get_nonhwe_snps(filename)
    print(nonhwe_snps)
