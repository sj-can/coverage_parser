from __future__ import print_function
import os, subprocess
from program import Argument_handler
from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises

def setup_module():
    print(__name__, ': setup_module() ~~~~~~~~~~~~~~~~~~~~~~')
 
def teardown_module():
    print(__name__, ': teardown_module() ~~~~~~~~~~~~~~~~~~~')

#creates the TestArgument_handler object, sets up the class method then set up for each test, then the test, then the teardowns. 
class TestArgument_handler(object):

    @classmethod
    def setup_class(cls):
	print(__name__, ': TestArgument_handler.setup_class() ----------')


    @classmethod
    def teardown_class(cls):
	print(__name__, ': TestArgument_handler.teardown_class() -------')
	#remove intemediate files
        for file in os.listdir(os.getcwd()):
            if file.startswith("NM_") or file.endswith("plottable.coverage"):
                remove_command = ["rm", file]
                subprocess.call(remove_command)
	os.remove('WE0345_KCNJ11.png')

    #occurs before each test method
    def setUp(self):
	with open('test_exome_list', 'w') as test_exomes:
	    test_exomes.write('WE0345')
	with open('test_gene_list', 'w') as test_gene:
	    test_gene.write('KCNJ11')

    #occurs after each test method
    def teardown(self):
	os.remove('test_exome_list')
	os.remove('test_gene_list')

    def test_init(self):
	A = Argument_handler('-s', 'test_exome_list', '-g', 'test_gene_list')
        assert_equal(A.arg_dict, {'sample_names': 'test_exome_list', 'gene_names': 'test_gene_list'})

    def test_create_parser(self):
	A = Argument_handler()
        this_parser = A.create_parser()
	parsed_args = this_parser.parse_args(['-s', 'test_exome_list', '-g', 'test_gene_list'])
	print(parsed_args)
        assert_equal(parsed_args.sample_names, 'test_exome_list')
	assert_equal(parsed_args.gene_names, 'test_gene_list')

    #end_to_end_test
    def test_handler(self):
        A =  Argument_handler('-s', 'test_exome_list', '-g', 'test_gene_list')
        A.handler()
        expected_transcript = os.path.isfile('NM_000525.3_reverse')
	assert_equal(expected_transcript, True)
	if expected_transcript:
	    command = ['wc', '-l', 'WE0345_KCNJ11.plottable.coverage']
	    process = subprocess.Popen(command, stdout=subprocess.PIPE)
	    line_count = process.communicate()[0]
	    assert_equal(line_count, '2957 WE0345_KCNJ11.plottable.coverage\n')
	    command_2 = ['wc', '-l', 'NM_001166290.1.intervals.exons']
	    process_2 = subprocess.Popen(command_2, stdout=subprocess.PIPE)
	    line_count_2 = process_2.communicate()[0]
	    assert_equal(line_count_2, '2957 NM_001166290.1.intervals.exons\n')
	   
	    


	
	

