
## CLASSES 

class RHTranscript(object):
    def __init__(self, line = None):
        if line:
            self.parse_from_line(line)
        else:
            self.txnstart = None
            self.tlnstart = None
            self.tlnend = None
            self.txnend = None
            self.strand = None
            self.name = None
        self.operons = []

    def __len__(self):
        return self.txnend - self.txnstart

    def parse_from_line(self, line):
        linearr = line.rstrip().split("\t")
        self.strand = linearr[4]
        if self.strand == "-":
            self.txnend = type_default(linearr[0], int)
            self.tlnend = type_default(linearr[1], int)
            self.tlnstart = type_default(linearr[2], int)
            self.txnstart = type_default(linearr[3], int)
        else:
            if self.strand != "+":
                self.strand = "."
            self.txnstart = type_default(linearr[0], int)
            self.tlnstart = type_default(linearr[1], int)
            self.tlnend = type_default(linearr[2], int)
            self.txnend = type_default(linearr[3], int)
        self.name = linearr[5]

    def add_operon(self, operon):
        self.operons.append(operon)

    def to_bed(self, chrm):
        start = self.txnstart if self.txnstart else self.tlnstart 
        end = self.txnend if self.txnend else self.tlnend
        # rockhopper is in 1-based coordinates. Must subtract 1 from the
        # start when converting to bed
        return "%s\t%d\t%d\t%s\t.\t%s"%(chrm, start - 1, end, self.name,self.strand)
        


class RHTranscripts(object):
    def __init__(self, infile = None):
        self.data = []
        if infile:
            self.parse_from_file(infile)

    def __iter__(self):
        for datum in self.data:
            yield datum

    def parse_from_file(self, infile):
        with open(infile) as inf:
            # skip the header
            inf.readline()
            for line in inf:
                if len(line) >= 5:
                    self.data.append(RHTranscript(line))

    def add_transcript(self, transcript):
        self.data.append(transcript)

class RHOperon(object):
    def __init__(self, line = None):
        if line:
            self.parse_from_line(line)
        else:
            self.start = None
            self.end = None
            self.strand = None
            self.genes = None
        self.transcripts = []

    def parse_from_line(self, line):
        linearr = line.rstrip().split("\t")
        self.start = int(linearr[0])
        self.end = int(linearr[1])
        self.strand = linearr[2]
        self.genes = linearr[4].rstrip().split(", ")

    def add_transcript(self, transcript):
        self.transcripts.append(transcript)

    def connect_transcripts(self, transcripts):
        for transcript in transcripts:
            if transcript.name in self.genes:
                self.add_transcript(transcript)
                transcript.add_operon(self)

    def update_by_transcripts(self):
        if self.transcripts:
            all_starts = [txn.txnstart for txn in self.transcripts]
            all_starts.append(self.start)
            all_ends = [txn.txnend for txn in self.transcripts]
            all_ends.append(self.end)
            self.start = min(filter(lambda x: x is not None, all_starts))
            self.end = max(filter(lambda x: x is not None, all_ends))

    def to_transcript(self):
        new_transcript = RHTranscript()
        new_transcript.txnstart = self.start
        new_transcript.txnend = self.end
        new_transcript.strand = self.strand
        new_transcript.name = "&".join(self.genes)
        return new_transcript

class RHOperons(object):
    def __init__(self, infile = None):
        self.data = []
        if infile:
            self.parse_from_file(infile)

    def __iter__(self):
        for datum in self.data:
            yield datum

    def parse_from_file(self, infile):
        with open(infile) as inf:
            # skip the header
            inf.readline()
            for line in inf:
                self.data.append(RHOperon(line))

    def match_to_transcripts(self, transcripts):
        for operon in self:
            operon.connect_transcripts(transcripts)

    def update_by_transcripts(self):
        for operon in self:
            if operon.transcripts:
                operon.update_by_transcripts()

## FUNCTIONS

def type_default(val, type_func, default_val = None):
    try:
        outval = type_func(val)
    except ValueError:
        outval = default_val
    return outval

def consolidate_transcripts_and_operons(operons, transcripts, chrm_name, min_size = 30):
    # match operons to transcripts and vice-versa
    operons.match_to_transcripts(transcripts)
    operons.update_by_transcripts()
    # get all the operons as transcripts
    final_transcripts = RHTranscripts()
    for operon in operons:
        final_transcripts.add_transcript(operon.to_transcript())
    for transcript in transcripts:
        # if already accounted for in the operons then don't add it
        if len(transcript.operons) > 0:
            continue
        # filter out any newly found RNAs less than 30 bp as they are 
        # probably not anything.
        if transcript.name == "-" and len(transcript) < min_size:
            continue
        else:
            final_transcripts.add_transcript(transcript)
    # update any duplicate names
    unique_names = {}
    for transcript in final_transcripts:
        # make a dictionary of all unique names
        if transcript.name == "-":
            transcript.name = chrm_name + "_predicted_RNA_%s"%(transcript.strand)
        this_list = unique_names.get(transcript.name, [])
        this_list.append(transcript)
        unique_names[transcript.name] = this_list
    # update names to be unique
    for entry in unique_names:
        these_transcripts = unique_names[entry]
        if len(these_transcripts) > 1:
            for i, transcript in enumerate(these_transcripts):
                transcript.name = transcript.name + ";%i"%(i)
    return final_transcripts



if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser(description="Consolidate transcripts and operon output from Rockhopper into a unique set of transcripts")
    parser.add_argument("operons_file", type=str, help="operons file input")
    parser.add_argument("transcripts_file", type=str, help="transcripts file input")
    parser.add_argument("chrm_name", type=str, help="Name of chromosome")
    parser.add_argument("--min_size", type=int, default=30, help="minimum length allowed for a newly discovered RNA")
    parser.add_argument("--skip_unstranded", action = "store_true", help="skip unstranded features?")
    args = parser.parse_args()

    operons = args.operons_file
    transcripts = args.transcripts_file
    chrm_name = args.chrm_name

    all_operons = RHOperons(operons)
    all_transcripts = RHTranscripts(transcripts)
    consolidated = consolidate_transcripts_and_operons(all_operons, all_transcripts, chrm_name, min_size = args.min_size)
    if args.skip_unstranded:
        for transcript in consolidated:
            if transcript.strand != ".": 
                sys.stdout.write(transcript.to_bed(chrm_name) + "\n")
    else:
        for transcript in consolidated:
            sys.stdout.write(transcript.to_bed(chrm_name) + "\n")

