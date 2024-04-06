import click
from vpolo.alevin import parser

@click.command()
@click.option('--one',  help='path to first alevin')
@click.option('--two',  help='path to second alevin')
def match(one, two):
	orig = parser.read_quants_bin(one)
	new = parser.read_quants_bin(two)
	diff = orig - new
	diff_sum = diff.sum().sum()
	if diff_sum > 0:
		print("Test failed with diff: {}".format(diff_sum))
	else:
		print("Test Passed")

	orig_new = parser.read_quants_bin(one, clipped='True')
	new_new = parser.read_quants_bin(two, clipped='True')
	print("Number of zerod cell in one: {}".format( len(orig.index)-len(orig_new.index) ))
	print("Number of zerod cell in two: {}".format( len(new.index)-len(new_new.index) ))


if __name__=="__main__":
    match()
