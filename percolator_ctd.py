import pprint
import os
import tempfile

from CTDopts.CTDopts import CTDModel, args_from_file, parse_cl_directives, flatten_dict, override_args, ArgumentRestrictionError

#The tool's model
msgf_percolator_model = CTDModel(
    name='Start_Percolator',
    version='0.1',
    description='This tool converts a mzid file to a tab file, starts the percolator and converts the percolators pout '
                'file to mzid')

msgf_percolator_model.add(
    name='ctd',
    type=str,
    default='false',
    description='Path to CTD file if already available'
)

######## mzid2pin ################

mzid2pinparams = msgf_percolator_model.add_group('mzid2pin', 'Grouped settings')

mzid2pinparams.add(
    'mzid2pin_path',
    required=True,
    type=str,
    description='set the path to your mzid2pin binary'
)

mzid2pinparams.add(
    'outputTab',
    required=True,
    type=str,
    default='--outputTab pin.tab',
    description='save output in a tab delimited file'
)

mzid2pinparams.add(
    'outputXML',
    type=str,
    description='save output in a pin-xml format file'
)

mzid2pinparams.add(
    'matches',
    type=int,
    description='Maximal number of matches to take in consideration'
)

mzid2pinparams.add(
    'verbose',
    type=int,
    num_range=(0, 5),
    default=2,
    description='Set verbosity of output: 0=no processing info, 5=all, default is: 2'
)

mzid2pinparams.add(
    'aa_freq',
    type=bool,
    description='This is a flag: Type 1 to set flag. calculate amino acid frequency features.'
)

mzid2pinparams.add(
    'PTM',
    type=bool,
    description='This is a flag: Type 1 to set flag. calculate feature for number of post-translational modifications'
)

mzid2pinparams.add(
    'enzyme',
    type=str,
    choices=['no_enzyme', 'elastase', 'pepsin', 'proteinasek', 'thermolysin', 'chymotrypsin', 'lys-n', 'lys-c', 'arg-c', 'asp-n', 'glu-c', 'trypsin'],
    default='trypsin',
    description='type of enzyme'
)

mzid2pinparams.add(
    'PNGaseF',
    type=bool,
    description='This is a flag: Type 1 to set flag. Calculate feature based on N-linked glycosylation pattern resulting from a PNGaseF treatment. (N[*].[ST])'
)

mzid2pinparams.add(
    'ms2_file',
    type=str,
    description='File containing spectra and retention time. The file could be in mzXML, MS2 or compressed MS2 file.'
)

mzid2pinparams.add(
    'isotope',
    type=bool,
    description='This is a flag: Type 1 to set flag. Mass difference calculated to closest isotope mass rather than to the average mass.'
)

mzid2pinparams.add(
    'psm_annotation',
    type=str,
    description='An anotation scheme used to convert the psms from the search. An example if Q# was used to describe '
                'pyro-glu formation (UNIMOD:28), and S* and T* was used to describe phosphorylation (UNIMOD:21), '
                'we would use the option -p *:21:#:28'
)

mzid2pinparams.add(
    'pattern',
    type=str,
    description='Pattern used to identify the decoy PSMs'
)

mzid2pinparams.add(
    'databases',
    type='input-file',
    file_formats=['fasta'],
    description='Link to the fasta database/s used in the search against the spectra file/s <target.fasta,[decoy.fasta]>'
                ' (Including this option will add the proteins to the generated pin file).'
)

mzid2pinparams.add(
    'cleavages',
    type=int,
    num_range=(0, None),
    default=0,
    description='Number of allowed miss cleavages used in the search engine (default 0) '
                '(Only valid when using option --databases).'
)

mzid2pinparams.add(
    'min_length',
    type=int,
    num_range=(0, 20),
    default=6,
    description='Minimum peptide length allowed used in the search engine (default 6) '
                '(Only valid when using option --databases).'
)

mzid2pinparams.add(
    'max_length',
    type=int,
    num_range=(0, 100),
    default=40,
    description='Maximum peptide length allowed used in the search engine (default 40) '
                '(Only valid when using option --databases).'
)

mzid2pinparams.add(
    'min_mass',
    type=int,
    num_range=(0, None),
    default=400,
    description='Minimum peptide mass allowed used in the search engine (default 400) '
                '(Only valid when using option --databases).'
)

mzid2pinparams.add(
    'max_mass',
    type=int,
    num_range=(0, None),
    default=6000,
    description='Maximum peptide mass allowed used in the search engine (default 6000) '
                '(Only valid when using option --databases).'
)

mzid2pinparams.add(
    'target_input',
    required=True,
    type='input-file',
    file_formats=['mzid'],
    description='Target input file'
)

# the decoy database for the converter
mzid2pinparams.add(
    'decoy_input',
    required=True,
    type='input-file',
    file_formats=['mzid'],
    description='Decoy input file'
)

#### percolator ####
pout = tempfile.NamedTemporaryFile()
percolatorparams = msgf_percolator_model.add_group('percolator', 'Grouped settings')

percolatorparams.add(
    'percolator_path',
    required=True,
    type=str,
    description='set the path to your percolator binary'
)

percolatorparams.add(
    'xmloutput',
    type=str,
    required=True,
    default='-X pout.xml',
    description='Name of the output file pout. Default is already set to temp file pout. do not set this parameter!'
)

percolatorparams.add(
    'stdinput',
    type=bool,
    description='Read xml-input format (pin) from standard input'
)

percolatorparams.add(
    'decoy_xml_output',
    type=bool,
    description='Include decoys (PSMs, peptide and/or proteins) in the xml-output. Only available if -X is used '
                '(in our case always)'
)

percolatorparams.add(
    'Cpos',
    type=int,
    description='Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified'
)

percolatorparams.add(
    'Cneg',
    type=int,
    description='Cneg, penalty for mistakes made on positive examples. Set by cross validation if not specified '
                'or --Cpos not specified'
)

percolatorparams.add(
    'trainFDR',
    type=float,
    default=0.01,
    description='False discovery rate threshold to define positive examples in training. Set by cross validation if 0. '
                'Default is 0.01.'
)

percolatorparams.add(
    'testFDR',
    type=float,
    default=0.01,
    description='False discovery rate threshold for evaluating best cross validation result and the reported '
                'end result. Default is 0.01.'
)

percolatorparams.add(
    'maxiter',
    type=int,
    description='Maximal number of iterations'
)

percolatorparams.add(
    'quick_validation',
    type=float,
    default=0.01,
    description='False discovery rate threshold for evaluating best cross validation result and the reported '
                'end result. Default is 0.01.'
)

percolatorparams.add(
    'train_ratio',
    type=float,
    default=0.6,
    description='Fraction of the negative data set to be used as train set when only providing one negative set, '
                'remaining examples will be used as test set. Set to 0.6 by default.'
)

percolatorparams.add(
    'tab_out',
    type=str,
    description='Output the computed features to the given file in tab-delimited format. '
                'A file with the features with the given file name will be created'
)

percolatorparams.add(
    'xml_in',
    type=str,
    description='Input file given in the deprecated pin-xml format generated by e.g. sqt2pin with the -k option'
)

percolatorparams.add(
    'weights',
    type=str,
    description='Output final weights to the given file'
)

percolatorparams.add(
    'init_weights',
    type=str,
    description='Read initial weights from the given file (one per line)'
)

percolatorparams.add(
    'default_direction',
    type=str,
    description='The most informative feature given as the feature name, can be negated to indicate that a lower value is better.'
)

percolatorparams.add(
    'verbose',
    type=int,
    num_range=(0,5),
    default=2,
    description='Set verbosity of output: 0=no processing info, 5=all, default is 2'
)

percolatorparams.add(
    'unitnorm',
    type=bool,
    description='Use unit normalization [0-1] instead of standard deviation normalization'
)

percolatorparams.add(
    'test_each_iteration',
    type=bool,
    description='This is a flag: Type 1 to set flag. Measure performance on test set each iteration'
)

percolatorparams.add(
    'override',
    type=bool,
    description='This is a flag: Type 1 to set flag. Override error check and do not fall back on default score vector in case of suspect score vector'
)

percolatorparams.add(
    'seed',
    type=int,
    default=1,
    description='Setting seed of the random number generator. Default value is 1'
)

percolatorparams.add(
    'klammer',
    type=bool,
    description='This is a flag: Type 1 to set flag. Retention time features calculated as in Klammer et al.'
)

percolatorparams.add(
    'doc',
    type=bool,
    description='Include description of correct features.'
)

percolatorparams.add(
    'results',
    type=str,
    description='Output tab delimited results to a file instead of stdout'
)

percolatorparams.add(
    'decoy_results',
    type=str,
    description='Output tab delimited results for decoys into a file'
)

percolatorparams.add(
    'only_psms',
    type=bool,
    description='This is a flag: Type 1 to set flag. Output tab delimited results for decoys into a file'
)

percolatorparams.add(
    'no_schema_validation',
    type=bool,
    description='This is a flag: Type 1 to set flag. skip validation of input file against xml schema'
)

percolatorparams.add(
    'protein',
    type=bool,
    description='This is a flag: Type 1 to set flag. Output protein level probabilities'
)

percolatorparams.add(
    'fido_alpha',
    type=float,
    description='Probability with which a present protein emits an associated peptide '
                '(to be used jointly with the --protein option) Set by grid search if not specified.'
)

percolatorparams.add(
    'fido_beta',
    type=float,
    description='Probability of the creation of a peptide from noise '
                '(to be used jointly with the --peptide option). Set by grid search if not specified'
)

percolatorparams.add(
    'fido_gamma',
    type=float,
    description='Prior probability of that a protein is present in the sample '
                '(to be used with the --protein option). Set by grid search if not specified'
)

percolatorparams.add(
    'allow_protein-group',
    type=bool,
    description='This is a flag: Type 1 to set flag. Treat ties as if it were one protein (Only valid if option --protein is active).'
)

percolatorparams.add(
    'protein_level_pi0',
    type=bool,
    description='use pi_0 value when calculating empirical q-values '
                '(no effect if option --fido-protein-group-level-interference is activated) '
                '(Only valid if option --protein is active).'
)

percolatorparams.add(
    'empirical_protein_q',
    type=bool,
    description='This is a flag: Type 1 to set flag. output empirical q-values and p-values (from target-decoy analysis) '
                '(Only valid if option --protein is active).'
)

percolatorparams.add(
    'fido_no_group_proteins',
    type=bool,
    description='This is a flag: Type 1 to set flag. disactivates the grouping of proteins with similar connectivity, '
                'for example if proteins P1 and P2 have the same peptides matching both of them, '
                'P1 and P2 will not be grouped as one protein. '
                '(Only valid if option --protein is active).'
)

percolatorparams.add(
    'fido_no_separate_proteins',
    type=bool,
    description='This is a flag: Type 1 to set flag. Proteins graph will not be separated in sub-graphs '
                '(Only valid if option --protein is active).'
)

percolatorparams.add(
    'fido_no_prune_proteins',
    type=bool,
    description='This is a flag: Type 1 to set flag. it does not prune peptides with a verylow score (~0.0) which means that if a peptide '
                'with a very low score is matching two proteins, when we prune the peptide, '
                'it will be duplicated to generate two new protein groups '
                '(Only valid if option --protein is active).'
)

percolatorparams.add(
    'fido_gridsearch_depth',
    type=int,
    default=0,
    description='Setting depth 0 or 1 or 2 from low depth to high depth(less computational time) '
                'of the grid search for the estimation Alpha,Beta and Gamma parameters for fido '
                '(Only valid if option --protein is active). Default value is 0'
)

percolatorparams.add(
    'pattern',
    type=str,
    default='random',
    description='Define the text pattern to identify the decoy proteins and/or PSMs, '
                'set this up if the label that idenfifies the decoys in the database '
                'is not the default (by default : random) (Only valid if option --protein  is active).'
)

percolatorparams.add(
    'fido_reduce_tree_in_gridsearch',
    type=bool,
    description='This is a flag: Type 1 to set flag. Reduce the tree of proteins (removing low scored proteins) in order to estimate '
                'alpha,beta and gamma faster.(Only valid if option --protein is active).'
)

percolatorparams.add(
    'post_processing_tdcn',
    type=bool,
    description='This is a flag: Type 1 to set flag. Use target decoy competition to compute peptide probabilities. '
                '(recommended when using --protein).'
)

percolatorparams.add(
    'grid_search_mse_threshold',
    type=float,
    description='Q-value threshold that will be used in the computation of the MSE and ROC AUC score '
                'in the grid search (recommended 0.05 for normal size datasets and 0.1 for big size datasets).'
                '(Only valid if option --protein is active).'
)

percolatorparams.add(
    'fido_truncation',
    type=bool,
    description='This is a flag: Type 1 to set flag. Proteins with a very low score (< 0.001) will be truncated (assigned 0.0 probability). '
                '(Only valid if option --protein is active).'
)

percolatorparams.add(
    'fido_protein_group_level_inference',
    type=bool,
    description='This is a flag: Type 1 to set flag. Uses protein group level inference, each cluster of proteins is either present or not, '
                'therefore when grouping proteins discard all possible combinations for each group. '
                '(Only valid if option --protein is active and --fido-no-group-proteins is inactive).'
)

percolatorparams.add(
    'input_file', # the tab file for the percolator
    required=True,
    type=str,
    default='pin.tab',
    description='The percolator tool needs a .tab file, which is set per default. Do not set this parameter, '
                'because it is a tempfile!'
)


###### pout2mzid ######

pout2mzidparams = msgf_percolator_model.add_group('pout2mzid', 'Grouped settings')

pout2mzidparams.add(
    'pout2mzid_path',
    required=True,
    type=str,
    description='set the path to your pout2mzid binary'
)

pout2mzidparams.add(
    'percolatorfile',
    type='input-file',
    default='-p pout.xml',
    description='Percolator Out XML result file. Do net set this parameter, because a temporary file is used'
)

pout2mzidparams.add(
    'mzidfile',
    type='input-file',
    required=True,
    description='MzIdentML input file'
)

pout2mzidparams.add(
    'output',
    type='input-file',
    description='Outputs the results to original filename+FILE+.mzid'
)

pout2mzidparams.add(
    'filesmzid',
    type='input-file',
    description='File containing a list of mzIdentML filenames'
)

pout2mzidparams.add(
    'decoy',
    type=str,
    description='Only adds results to entries with decoy set to true. DEFAULT: false'
)

pout2mzidparams.add(
    'validate',
    type=bool,
    description='This is a flag: Type 1 to set flag. Sets that validation of XML schema should not be performed. Faster parsing.'
)

pout2mzidparams.add(
    'warning',
    type=bool,
    description='This is a flag: Type 1 to set flag. Sets that upon warning the software should terminate'
)

args = msgf_percolator_model.parse_cl_args()

try:
    args = args_from_file(args['ctd'])
    print "try to load ctd"

    ### MZID2PIN Parameter ###
    try:
        args['mzid2pin']['mzid2pin_path']
        args['mzid2pin']['target_input']
        args['mzid2pin']['decoy_input']
    except:
        print "Path or decoy/target-files missing!"
    str = ''
    try:
        if args['mzid2pin']['outputXML'] != '':
            str += ' --outputXML \'%s\'' % args['mzid2pin']['outputXML']
    except:
        pass
    try:
        if args['mzid2pin']['matches'] != '':
            str += ' --matches %s' % args['mzid2pin']['matches']
    except:
        pass
    try:
        if args['mzid2pin']['verbose'] != '2':
            str += ' --verbose %s' % args['mzid2pin']['verbose']
    except:
        pass
    try:
        if args['mzid2pin']['aa_freq'] == 'true':
            str += ' --aa-freq'
    except:
        pass
    try:
        if args['mzid2pin']['PTM'] == 'true':
            str += ' --PTM'
    except:
        pass
    try:
        if args['mzid2pin']['enzyme'] != 'trypsin':
            str += ' --enzyme \'%s\'' % args['mzid2pin']['enzyme']
    except:
        pass
    try:
        if args['mzid2pin']['PNGaseF'] == 'true':
            str += ' --PNGaseF'
    except:
        pass
    try:
        if args['mzid2pin']['ms2_file'] != '':
            str += ' --ms2-file \'%s\'' % args['mzid2pin']['ms2_file']
    except:
        pass
    try:
        if args['mzid2pin']['isotope'] == 'true':
            str += ' --isotope'
    except:
        pass
    try:
        if args['mzid2pin']['psm_annotation'] != '':
            str += ' --psm-annotation \'%s\'' % args['mzid2pin']['psm_annotation']
    except:
        pass
    try:
        if args['mzid2pin']['pattern'] != '':
            str += ' --pattern \'%s\'' % args['mzid2pin']['pattern']
    except:
        pass
    try:
        if args['mzid2pin']['databases'] != '':
            str += ' --databases \'%s\'' % args['mzid2pin']['databases']
    except:
        pass
    try:
        if args['mzid2pin']['cleavages'] != '0':
            str += ' --cleavages %s' % args['mzid2pin']['cleavages']
    except:
        pass
    try:
        if args['mzid2pin']['min_length'] != '6':
            str += ' --min-length %s' % args['mzid2pin']['min_length']
    except:
        pass
    try:
        if args['mzid2pin']['max_length'] != '40':
            str += ' --max-length %s' % args['mzid2pin']['max_length']
    except:
        pass
    try:
        if args['mzid2pin']['min_mass'] != '400':
            str += ' --min-mass %s' % args['mzid2pin']['min_mass']
    except:
        pass
    try:
        if args['mzid2pin']['max_mass'] != '6000':
            str += ' --max-mass %s' % args['mzid2pin']['max_mass']
    except:
        pass
    syscall_mzid = args['mzid2pin']['mzid2pin_path'] + ' ' + args['mzid2pin']['outputTab'] + str + ' ' + args['mzid2pin']['target_input'] + ' ' + args['mzid2pin']['decoy_input']


    ### PERCOLATOR PARAMETER
    try:
        args['percolator']['percolator_path']
        args['percolator']['input_file']
        args['percolator']['xmloutput']
    except:
        print "Pfad, input oder output fehlen!"
    str = ''
    try:
        if args['percolator']['input_file'] != 'pin.tab':
            str += ' --input-file %s' % args['percolator']['input_file']
    except:
        pass
    try:
        if args['percolator']['stdinput'] == 'true':
            str += ' --stdinput'
    except:
        pass
    try:
        if args['percolator']['decoy_xml_output'] == 'true':
            str += ' --decoy-xml-output'
    except:
        pass
    try:
        if args['percolator']['Cpos'] != '':
            str += ' --Cpos %s' % args['percolator']['Cpos']
    except:
        pass
    try:
        if args['percolator']['Cneg'] != '':
            str += ' --Cneg %s' % args['percolator']['Cneg']
    except:
        pass
    try:
        if args['percolator']['trainFDR'] != '0.01':
            str += ' --trainFDR %s' % args['percolator']['trainFDR']
    except:
        pass
    try:
        if args['percolator']['testFDR'] != '0.01':
            str += ' --testFDR %s' % args['percolator']['testFDR']
    except:
        pass
    try:
        if args['percolator']['maxiter'] != '':
            str += ' --maxiter %s' % args['percolator']['maxiter']
    except:
        pass
    try:
        if args['percolator']['quick_validation'] == 'true':
            str += ' --quick-validation'
    except:
        pass
    try:
        if args['percolator']['train_ratio'] != '0.6':
            str += ' --train-ratio %s' % args['percolator']['train_ratio']
    except:
        pass
    try:
        if args['percolator']['tab_out'] != '':
            str += ' --tab-out \'%s\'' % args['percolator']['tab_out']
    except:
        pass
    try:
        if args['percolator']['xml_in'] != '':
            str += ' --xml-in \'%s\'' % args['percolator']['xml_in']
    except:
        pass
    try:
        if args['percolator']['weights'] != '':
            str += ' --weights \'%s\'' % args['percolator']['weights']
    except:
        pass
    try:
        if args['percolator']['init_weights'] != '':
            str += ' --init-weights \'%s\'' % args['percolator']['init_weights']
    except:
        pass
    try:
        if args['percolator']['default_direction'] != '':
            str += ' --default-direction \'%s\'' % args['percolator']['default_direction']
    except:
        pass
    try:
        if args['percolator']['verbose'] != '2':
            str += ' --verbose %s' % args['percolator']['verbose']
    except:
        pass
    try:
        if args['percolator']['unitnorm'] == 'true':
            str += ' --unitnorm'
    except:
        pass
    try:
        if args['percolator']['test_each_iteration'] == 'true':
            str += ' --test-each-iteration'
    except:
        pass
    try:
        if args['percolator']['override'] == 'true':
            str += ' --override'
    except:
        pass
    try:
        if args['percolator']['seed'] != '1':
            str += ' --seed %s' % args['percolator']['seed']
    except:
        pass
    try:
        if args['percolator']['klammer'] == 'true':
            str += ' --klammer'
    except:
        pass
    try:
        if args['percolator']['doc'] == 'true':
            str += ' --doc'
    except:
        pass
    try:
        if args['percolator']['results'] != '':
            str += ' --results %s' %args['percolator']['results']
    except:
        pass
    try:
        if args['percolator']['decoy_results'] != '':
            str += ' --decoy-results %s' %args['percolator']['decoy_results']
    except:
        pass
    try:
        if args['percolator']['only_psms'] == 'true':
            str += ' --only-psms'
    except:
        pass
    try:
        if args['percolator']['no_schema_validation'] == 'true':
            str += ' --no-schema-validation'
    except:
        pass
    try:
        if args['percolator']['protein'] == 'true':
            str += ' --protein'
    except:
        pass
    try:
        if args['percolator']['fido_alpha'] != '':
            str += ' --fido-alpha %s' %args['percolator']['fido_alpha']
    except:
        pass
    try:
        if args['percolator']['fido_beta'] != '':
            str += ' --fido-beta %s' %args['percolator']['fido_beta']
    except:
        pass
    try:
        if args['percolator']['fido_gamma'] != '':
            str += ' --fido-gamma %s' %args['percolator']['fido_gamma']
    except:
        pass
    try:
        if args['percolator']['allow_protein_group'] == 'true':
            str += ' --allow-protein-group'
    except:
        pass
    try:
        if args['percolator']['protein_level_pi0'] == 'true':
            str += ' --protein-level-pi0'
    except:
        pass
    try:
        if args['percolator']['empirical_protein_q'] == 'true':
            str += ' --empirical-protein-q'
    except:
        pass
    try:
        if args['percolator']['fido_no_group_proteins'] == 'true':
            str += ' --fido-no-group-proteins'
    except:
        pass
    try:
        if args['percolator']['fido_no_separate_proteins'] == 'true':
            str += ' --fido-no-Separate-proteins'
    except:
        pass
    try:
        if args['percolator']['fido_no_prune_proteins'] == 'true':
            str += ' --fido-no-prune-proteins'
    except:
        pass
    try:
        if args['percolator']['fido_gridsearch_depth'] == '':
            str += ' --fido-gridsearch-depth %s' %args['percolator']['fido_gridsearch_depth']
    except:
        pass
    try:
        if args['percolator']['pattern'] != 'random':
            str += ' --pattern %s' %args['percolator']['pattern']
    except:
        pass
    try:
        if args['percolator']['fido_reduce_tree_in_gridsearch'] == 'true':
            str += ' --fido-reduce-tree-in-gridsearch'
    except:
        pass
    try:
        if args['percolator']['post_processing_tdcn'] == 'true':
            str += ' --post-processing-tdcn'
    except:
        pass
    try:
        if args['percolator']['grid_search_mse_threshold'] == 'true':
            str += ' --grid-search-mse-threshold'
    except:
        pass
    try:
        if args['percolator']['fido_truncation'] == 'true':
            str += ' --fido-truncation'
    except:
        pass
    try:
        if args['percolator']['fido_protein_group_level_inference'] == 'true':
            str += ' --fido-protein-group-level-inference'
    except:
        pass
    syscall_perc = args['percolator']['percolator_path'] + ' ' + args['percolator']['xmloutput'] + str + ' ' + args['percolator']['input_file']

    ### POUT2MZID PARAMETER
    try:
        args['pout2mzid']['pout2mzid_path']
        args['pout2mzid']['percolatorfile']
        args['pout2mzid']['mzidfile']
    except:
        print "Path, pout or mzid missing!"
    str = ''
    try:
        if args['pout2mzid']['output'] != '':
            str += ' --output \'%s\'' % args['pout2mzid']['output']
    except:
        pass
    try:
        if args['pout2mzid']['filesmzid'] != '':
            str += ' --filesmzid \'%s\'' % args['pout2mzid']['filesmzid']
    except:
        pass
    try:
        if args['pout2mzid']['decoy'] == 'true':
            str += ' --decoy'
    except:
        pass
    try:
        if args['pout2mzid']['validate'] == 'true':
            str += ' --validate'
    except:
        pass
    try:
        if args['pout2mzid']['warning'] == 'true':
            str += ' --warning'
    except:
        pass
    syscall_pout = args['pout2mzid']['pout2mzid_path'] + ' ' + args['pout2mzid']['percolatorfile'] + ' -m ' + args['pout2mzid']['mzidfile'] + str
    os.system('%s' % syscall_mzid)
    os.system('%s' % syscall_perc)
    os.system('%s' % syscall_pout)
except:
    print 'successfully wrote ctd'
    cl_args = msgf_percolator_model.parse_cl_args()
    write_cl_params = override_args(cl_args)
    validate = msgf_percolator_model.validate_args(write_cl_params)
    msgf_percolator_model.write_ctd('msgfpercolator.ctd', validate)