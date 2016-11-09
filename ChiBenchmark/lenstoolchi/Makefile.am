#Main Makefile for Lenstool

AUTOMAKE_OPTIONS        = foreign 
SUBDIRS 		= liblt src table_src utils 
EXTRA_DIST		= doc examples perl include examples_table

TABLE_DIR=./table_src

table: 
	cd $(TABLE_DIR) ; rm -f ./lenstool.tab ; ./lenstool_tab default.in
		
