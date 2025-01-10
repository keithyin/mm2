import pysam
from tqdm import tqdm


def main():
    fname = "/data/ccs_data/ccs_eval2024q4/20240607_Sync_H84_01_H01_Run0002_adapter-metric/20240607_Sync_H84_01_H01_Run0002_adapter.aligned.bam"
    with pysam.AlignmentFile(fname, mode="rb") as fh:
        cnt = 0
        for record in tqdm(fh.fetch(until_eof=True), desc=f"reading {fname}"):

            if record.is_reverse:
                cnt += 1
                print(record.query_sequence)
                print()

            if cnt > 10:
                break
            # print(
            #     "qname:{}, refstart:{}, refend:{}, qstart:{}, qend:{}. secondary:{}".format(
            #         record.query_name,
            #         record.reference_start,
            #         record.reference_end,
            #         record.query_alignment_start,
            #         record.query_alignment_end,
            #         record.is_secondary,
            #     )
            # )

    pass


if __name__ == "__main__":
    main()
