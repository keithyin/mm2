import mappy as mp
import pysam

#        "/data/ccs_data/ccs_eval2024q3/jinpu/ref_Saureus_ATCC25923.m.new.corrected.fasta"


def main():
    a = mp.Aligner(
        "/data/ccs_data/ccs_eval2024q3/jinpu/ref_dup.fasta"
    )  # load or build index
    if not a:
        raise Exception("ERROR: failed to load/build index")
    with pysam.AlignmentFile(
        "/data/ccs_data/ccs_eval2024q3/jinpu/20240711_Sync_Y0006_02_H01_Run0001_called.bam",
        mode="rb",
        check_sq=False,
    ) as in_bam:
        for record in in_bam.fetch(until_eof=True):
            seq = record.query_sequence
            print(seq)
            print("seq.name", record.query_name, "seq_len", len(seq))
            for hit in a.map(seq):
                print(hit, hit.is_primary)
                # print(
                #     "{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str)
                # )

            break


if __name__ == "__main__":

    main()
