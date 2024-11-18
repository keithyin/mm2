import pysam


def main():
    fname = "aligned.bam"
    with pysam.AlignmentFile(fname, mode="rb") as fh:
        for record in fh.fetch(until_eof=True):
            print(
                "qname:{}, refstart:{}, refend:{}, qstart:{}, qend:{}. secondary:{}".format(
                    record.query_name,
                    record.reference_start,
                    record.reference_end,
                    record.query_alignment_start,
                    record.query_alignment_end,
                    record.is_secondary,
                )
            )

    pass


if __name__ == "__main__":
    main()
