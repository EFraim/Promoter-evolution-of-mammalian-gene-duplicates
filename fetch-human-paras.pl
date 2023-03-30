use 5.12.0;
use warnings;

use Bio::EnsEMBL::Registry;
use Data::Dumper;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -db_version => 98
    );

my $c = 0;

##Adaptors
my $tree_adaptor        = $registry->get_adaptor("Multi", "compara", "GeneTree");
my $gene_member_adaptor = $registry->get_adaptor("Multi", "compara", "GeneMember");
my $caffe_adaptor = $registry->get_adaptor("Multi", "compara", "CAFEGeneFamily");
my $genome_db_adaptor   = Bio::EnsEMBL::Registry->get_adaptor( "Multi", "compara", "GenomeDB");

##Genome DBs
my @genomes = ('homo_sapiens'); #'anolis_carolinensis');

my @genome_dbs;
my $genome = shift @genomes;
my $genome_db = $genome_db_adaptor->fetch_by_name_assembly($genome);
my $genes = $gene_member_adaptor->fetch_all_by_GenomeDB($genome_db);
my %tree_per_gene = ();
@$genes  = @$genes[15400 .. $#$genes];
foreach my $gene (@$genes) {
    next if ($gene->biotype_group() ne 'coding');
    my $tree = $tree_adaptor->fetch_default_for_Member($gene);
    next unless $tree;
    print($tree->nhx_format("gene_id"));
    print "\n";
}
