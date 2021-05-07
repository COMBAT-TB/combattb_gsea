#!/usr/bin/env python
# Copyright 2021 SANBI, University of the Western Cape
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

from operator import itemgetter
import sys
import click
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests

from neo4j import GraphDatabase


class CombatTbGSEA:
    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self.driver.close()

    @staticmethod
    def _count_genes_tx(tx):
        result = tx.run("MATCH (g:Gene) return count(g)")
        return result.single()[0]

    def count_genes(self):
        with self.driver.session() as session:
            count = session.read_transaction(self._count_genes_tx)
        return count

    @staticmethod
    def _get_go_table_tx(tx, filter=""):
        query = (
            "MATCH (t:GOTerm) <-[:ASSOCIATED_WITH]- "
            "(:Protein) <-[:ENCODES]- (g:Gene) "
            + filter
            + " RETURN t.accession, t.name, count(g) "
            "AS t_count ORDER BY t_count DESC"
        )
        result = tx.run(query)
        return result.values()

    def get_goterm_prevalence(self, gene_list=None):
        if gene_list:
            filter = (
                "WHERE g.uniquename IN ["
                + ",".join([f'"{gene}"' for gene in gene_list])
                + "]"
            )
        else:
            filter = ""
        with self.driver.session() as session:
            values = session.read_transaction(self._get_go_table_tx, filter)
        return values


def enrichment_analysis(geneset, mode="over", multipletest_method="fdr_bh"):
    if mode not in ("over", "under"):
        raise ValueError("mode must be 'over' or 'under'")
    if multipletest_method not in ("fdr_bh", "bonferroni"):
        raise ValueError(
            "multipletest_method must be 'fdr_bh' (Benjamini-Hochberg) or 'bonferroni' (Bonferroni)"
        )
    db = CombatTbGSEA("bolt://neodb.sanbi.ac.za:7687", "", "")

    background_gene_count = db.count_genes()
    goterm_prevalence = db.get_goterm_prevalence()
    goterm_prevalence_dict = dict()
    for goterm_id, _, count in goterm_prevalence:
        goterm_prevalence_dict[goterm_id] = count

    gs_goterm_prevalence = db.get_goterm_prevalence(gene_list=geneset)
    geneset_size = len(geneset)
    p_vals = []
    if mode == "over":
        fe_test_mode = "greater"
    else:
        fe_test_mode = "less"
    for goterm_id, goterm_name, gs_goterm_count in gs_goterm_prevalence:
        # Fisher Exact test matrix
        #          GO term XXX      Not GO Term XXX
        # geneset     A                 B
        # background  C                 D
        background_goterm_count = goterm_prevalence_dict[goterm_id]
        fe_matrix = (
            (gs_goterm_count, geneset_size - gs_goterm_count),
            (background_goterm_count, background_gene_count - background_goterm_count),
        )
        _, p_val = fisher_exact(fe_matrix, alternative=fe_test_mode)
        p_vals.append((goterm_id, goterm_name, p_val))
    p_vals = sorted(p_vals, key=itemgetter(2))
    goterm_ids, goterm_names, uncorr_p_vals = zip(*p_vals)
    _, corrected_p_vals, _, _ = multipletests(
        uncorr_p_vals, is_sorted=True, method=multipletest_method
    )
    p_vals = zip(goterm_ids, goterm_names, uncorr_p_vals, corrected_p_vals)
    # for goterm_it, goterm_name, p_val, corr_p_val in p_vals:
    #     print(goterm_it, goterm_name, p_val, corr_p_val)
    return p_vals


@click.command()
@click.option("--mode", "-M", default="over", type=click.Choice(("over", "under")))
@click.option(
    "--multipletest_corr", "-T", default="BH", type=click.Choice(("BH", "BON"))
)
@click.argument("geneset_file", type=click.File("r"))
@click.argument("output_file", required=False, default=sys.stdout, type=click.File("w"))
def analyse_geneset(geneset_file, output_file, mode="over", multipletest_corr="BH"):
    geneset = []
    for line in geneset_file:
        geneset.append(line.rstrip().replace("rv_", "Rv"))
    if multipletest_corr == "BH":
        multi = "fdr_bh"
    else:
        multi = "bonferroni"
    results = enrichment_analysis(geneset, mode=mode, multipletest_method=multi)
    print(
        "GO Term ID\tGO Term Name\tCorrected p-val\tUncorrected p-val", file=output_file
    )
    for goterm_id, goterm_name, uncorr_p_val, corr_p_val in results:
        print(
            goterm_id, goterm_name, corr_p_val, uncorr_p_val, file=output_file, sep="\t"
        )
    output_file.close()


if __name__ == "__main__":
    analyse_geneset()
