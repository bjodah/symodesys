import argparse


class IntegrationArgumentParser(object):
    pass

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--at-src', action = 'store_true',
                    help = 'Create xyz in directory of source instead of at source')
parser.add_argument('-s', '--suffix', type = str, default = '',
                    help = 'Suffix for xyz file, e.g. `_final`')
parser.add_argument('-d', '--distort_normal', type = float, default = 0.0,
                    help = 'Distort along normal coordinate')
parser.add_argument('-n', '--normal_idx', type = int, default = -1,
                    help = 'Index of normal coordinate to distort along')
parser.add_argument('-g', '--geom_idx', type = int, default = -1,
                    help = 'Index of geometry to extract. Default: -1 (last)')

parser.add_argument('job_files',nargs='+',type=str,help='com- or log-files')
