use Bio::EnsEMBL::Registry;
use Bio::AlignIO;
use warnings;
use strict;
use Data::Dumper qw(Dumper);

my $input_file = $ARGV[0];
my $query_species = $ARGV[1];
my $other_species = $ARGV[2];
my $output_file = $ARGV[3];
my $version = $ARGV[4];

open(my $out_fh, '>', $output_file) or die "Could not open file '$output_file' $!";

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all("reg_conf.pl");

my $genome_db_adaptor = $registry->get_adaptor('Multi', 'compara', 'GenomeDB');
    
throw("Cannot connect to Compara") if (!$genome_db_adaptor);

my $genome_dbs;

foreach my $this_species ($query_species, $other_species) {
    my $this_meta_container_adaptor = $registry->get_adaptor(
        $this_species, 'core', 'MetaContainer');

    throw("Registry configuration file has no data for connecting to <$this_species>")
        if (!$this_meta_container_adaptor);

    my $this_production_name = $this_meta_container_adaptor->get_production_name;

    # Fetch Bio::EnsEMBL::Compara::GenomeDB object
    my $genome_db = $genome_db_adaptor->fetch_by_name_assembly($this_production_name);

    # Add Bio::EnsEMBL::Compara::GenomeDB object to the list
    push(@$genome_dbs, $genome_db);
}

my $slice_adaptor = $registry->get_adaptor(
    $query_species, 'core', 'Slice');
    
my $method_link_species_set_adaptor = $registry->get_adaptor(
    'Multi', 'compara', 'MethodLinkSpeciesSet');
    
my $mlss = $method_link_species_set_adaptor->fetch_by_method_link_type_GenomeDBs("LASTZ_NET", $genome_dbs);
    
my $genomic_align_block_adaptor = $registry->get_adaptor(
    'Multi', 'compara', 'GenomicAlignBlock');
    
my $alignIO = Bio::AlignIO->newFh(
    -interleaved => 0,
    -fh => $out_fh,
    -format => 'phylip',
    -idlength => 30
);

open(my $in_fh, '<', $input_file) or die "Cannot open $input_file: $!";

my $counter = 0;
while ( my $line = <$in_fh> ) {
	$counter = $counter + 1;
	print $counter, "\n";
    my @row = split( /\|/, $line );
    my $chrom = $row[0];
    my $start = $row[2];
	my $end = $row[3];

	print $out_fh '%', $line;
	
	my $query_slice = $slice_adaptor->fetch_by_region(
    	'toplevel',
    	$chrom,
    	$start,
    	$end);

	my $genomic_align_blocks =
    	$genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice(
      	$mlss,
      	$query_slice);
      
	my $all_aligns;

	foreach my $genomic_align_block( @{ $genomic_align_blocks }) {
	    my $restricted_gab = $genomic_align_block->restrict_between_reference_positions($start, $end);
		my $simple_align = $restricted_gab->get_SimpleAlign;
		push(@$all_aligns, $simple_align);
	}      
 
	foreach my $this_align (@$all_aligns) {
    	print $alignIO $this_align;
    	print $out_fh "|||\n";
	}
	
	print $out_fh "***\n";
}

close($in_fh);
close $out_fh;
