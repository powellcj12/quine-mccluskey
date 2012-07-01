quine-mccluskey
===============

## Introduction

This is a software tool which generates the minimized 2-level SOP (sum of products) and POS (product of sums) expressions of all single-output Boolean functions up to 10 literals.

## Usage

	make
	./QM

## Execution

First, you must specify how you will enter Boolean functions. They may be entered through the command line when running the progrm or in an input file that must be given.

Each line of input corresponds to 1 Boolean function. Each Boolean function must have at least one minterm and zero or more "don't cares". Below are some input examples:


	m(1,2,3,9,10)+(d5,7)
	m(0,2,3,5,6,7,8,10,11,14,15)
	m(1,5,3)+d(2,4)

Note that minterms do not need to be specified in sorted order.

## Output

The SOP and POS solutions will be printed to stdout.
