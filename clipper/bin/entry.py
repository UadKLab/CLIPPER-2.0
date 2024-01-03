import time
import logging

from typing import Dict, Optional, List, Tuple

from Bio import ExPASy, SwissProt


class Entry:
    def __init__(self, acc: str, seq: str) -> None:
        """
        Initializes an Entry object with a Uniprot accession and peptide sequence.

        Args:
            acc (str): Uniprot accession.
            seq (str): Peptide sequence.
        """
        self.acc = acc
        self.seq = seq
        self.record: Optional[SwissProt.Record] = None
        self.annot: Dict[str, str] = {}
        self.missing_uniprot_records = []
        self.warnings = {'retrieval': [],
                         'peptides not found': []}

    def get_record(self, sleep_time: float = 0.5) -> None:
        time.sleep(sleep_time)
        self.record = self._retrieve_uniprot_record()

    def _retrieve_uniprot_record(self) -> Optional[SwissProt.Record]:
        tries = 0
        while tries < 5:
            try:
                handle = ExPASy.get_sprot_raw(self.acc,)
                record = SwissProt.read(handle)
                return record
            except Exception as e:
                if tries == 4:
                    tries += 1
                    #logging.warning(f"Error retrieving Uniprot record for {self.acc}: {str(e)}")
                    self.warnings['retrieval'].append(f"Error retrieving Uniprot record for {self.acc}: {str(e)}")
                    return None
                else:
                    tries += 1
                    continue

    def parse_general(self) -> None:
        if not self.record:
            return

        self.annot["query_sequence"] = self.seq
        self.annot["query_accession"] = self.acc
        self.annot["name"] = self.record.entry_name
        self.annot["full_sequence"] = self.record.sequence
        self.annot["description"] = self.record.description
        self.annot["keywords"] = "|".join(self.record.keywords)

        codes, names = self._parse_go_annotations()
        if codes:
            self.annot["go_codes"] = "|".join(codes)
            self.annot["go_names"] = "|".join(names)

    def _parse_go_annotations(self) -> Tuple[List[str], List[str]]:
        codes = [ref[1] for ref in self.record.cross_references if ref[0] == "GO"]
        names = [ref[2] for ref in self.record.cross_references if ref[0] == "GO"]
        return codes, names

    def parse_cleavage(self, cleavagesitesize) -> None:
        if not self.record:
            return

        full_sequence = self.annot["full_sequence"]
        match = full_sequence.find(self.seq)

        if match == -1:
            self._handle_peptide_not_found(cleavagesitesize)
        else:
            self._handle_peptide_found(match, full_sequence, cleavagesitesize)

    def _handle_peptide_not_found(self, cleavagesitesize) -> None:
        
        self.annot["nterm_annot"] = "Not found"
        self.annot["start_pep"] = "Not found"
        self.annot["end_pep"] = "Not found"
        self.annot["cleavage_site"] = "Not found"
        self.annot[f"p{cleavagesitesize}_p{cleavagesitesize}prime"] = "Not found"
        self.cleavage_site = None
        self.warnings['peptides not found'].append(f"Peptide {self.seq} not found in protein sequence for {self.acc}.")

    def _handle_peptide_found(self, match: int, full_sequence: str, cleavagesitesize: int) -> None:
        
        pos = match
        self.cleavage_site = pos
        self.annot["p1_position"] = pos

        # generate cleavage site p4-p4', pad if not enough amino acids before cleavage
        if pos > cleavagesitesize - 1:
            preseq = full_sequence[pos - cleavagesitesize : pos]
            self.annot[f"p{cleavagesitesize}_p{cleavagesitesize}prime"] = preseq + self.seq[:cleavagesitesize]
        else:
            preseq = full_sequence[:pos]
            padding = "-" * (cleavagesitesize - pos)
            self.annot[f"p{cleavagesitesize}_p{cleavagesitesize}prime"] = padding + preseq + self.seq[:cleavagesitesize]

        self.preseq = full_sequence[:pos]

        self.annot["start_pep"] = pos + 1
        self.annot["end_pep"] = pos + len(self.seq)
        self.annot["cleavage_site"] = f"{preseq}({pos}).({pos+1}){self.seq}"
        self.annot["acc_length"] = len(full_sequence)

        # for now, only one signal peptide and propeptide are considered
        # if the proteins has more than one propeptides, the second is ignored
        sigpos = None
        propos = None
        transpos = None
        flag_sig = False
        flag_pro = False

        # find whether protein contains signal and propeptide sequences
        for feat in self.record.features:

            try:
                # last item of location list is up to last amino acid of e.g. signal peptide
                # so, if pos of peptide starts right after that, we need to increment by one
                feature_end = list(feat.location)[-1] + 1
                feature_start = list(feat.location)[0]

                if feat.type == "SIGNAL" and not flag_sig:
                    sigstart = feature_start
                    sigpos = feature_end
                    self.sigpos = sigpos
                    flag_sig = True
                elif feat.type == "PROPEP" and not flag_pro:
                    prostart = feature_start
                    propos = feature_end
                    self.propos = propos
                    flag_pro = True
                elif feat.type == "TRANSIT":
                    transtart = feature_start
                    transpos = feature_end
                    self.transpos = transpos

            except:
                feature_end = None

        # check where cleavage is located in the protein sequence
        # natural N-terminus or Met removal, no cleavage annotation necessary
        if pos == 0:
            self.annot["nterm_annot"] = "Met intact"
        elif pos == 1:
            self.annot["nterm_annot"] = "Met removed"

        # if both signal peptide and propeptide exist
        elif sigpos and propos:

            # normal cases in which 99% of proteins fall
            if pos == sigpos:
                self.annot["nterm_annot"] = "Signal removed"
            elif pos == propos:
                self.annot["nterm_annot"] = "Propeptide removed"
            elif sigstart < pos < sigpos:
                self.annot["nterm_annot"] = "Cleavage within signal peptide range"
            elif prostart < pos < propos:
                self.annot["nterm_annot"] = "Cleavage within propeptide range"
            elif pos > propos:
                self.annot["nterm_annot"] = "Internal"

            # special cases
            # if cleavage is before signal peptide, which is weird, annotate special case
            elif pos < sigstart:
                self.annot[
                    "nterm_annot"
                ] = "Internal, cleavage before signal peptide"

            # another weird case where cleavage is at the start of the signal peptide
            # assuming this is only in case the signal peptide is at C-term for some reason
            elif pos == sigstart:
                self.annot[
                    "nterm_annot"
                ] = "Signal removed, signal peptide at C-term or cleavage at start of singal peptide"

            # if propeptide is later in the sequence, cleavage is internal
            # if propeptide starts right after signal peptide, sigpos is equal to prostart
            # so this will not be executed as we covered that case on condition 3
            elif pos < prostart:
                self.annot["nterm_annot"] = "Internal"

            # in the case of sigpos not equal to sigstart, propeptide is in C-term
            elif pos == prostart:
                self.annot["nterm_annot"] = "Propeptide removed"

            # catch anything that could have slipped through, although we should have covered all cases
            else:
                self.annot["nterm_annot"] = "Exception"

        # if only signal peptide exists
        elif sigpos:

            if pos == sigpos:
                self.annot["nterm_annot"] = "Signal removed"
            elif sigstart < pos < sigpos:
                self.annot["nterm_annot"] = "Cleavage within signal peptide range"
            elif pos > sigpos:
                self.annot["nterm_annot"] = "Internal"

            # special cases
            elif pos < sigstart:
                self.annot[
                    "nterm_annot"
                ] = "Internal, cleavage before signal peptide"
            elif pos == sigstart:
                self.annot[
                    "nterm_annot"
                ] = "Signal removed, signal peptide at C-term or cleavage at start of singal peptide"
            else:
                self.annot["nterm_annot"] = "Exception"

        # if only propeptide exists
        elif propos:

            if pos == propos:
                self.annot["nterm_annot"] = "Propeptide removed"
            elif prostart < pos < propos:
                self.annot["nterm_annot"] = "Cleavage within propeptide range"
            elif pos > propos:
                self.annot["nterm_annot"] = "Internal"

            # special cases
            elif pos < prostart:
                self.annot["nterm_annot"] = "Internal"
            # propeptide is in the C-term
            elif pos == prostart:
                self.annot["nterm_annot"] = "Propeptide removed, C-term"
            else:
                self.annot["nterm_annot"] = "Exception"

        # if transit peptide exists, assumes transit peptides do not coexist with signal and propeptides
        elif transpos:

            if pos == transpos:
                self.annot["nterm_annot"] = "Transit peptide removed"
            elif transtart < pos < transpos:
                self.annot["nterm_annot"] = "Cleavage within transit peptide range"
            elif pos > transpos:
                self.annot["nterm_annot"] = "Internal"

        # if no signal or propeptide exist
        else:
            self.annot["nterm_annot"] = "Internal"

    def parse_protease(self):
        """Parses uniprot record to cheeck for PTM cleavage of substrate in the
        specified position"""

        proteases = []

        for feat in self.record.features:

            if feat.type == "SITE":
                # cleavage_site is p1' 0-index
                # location list is also 0-index [p1,p1']
                if list(feat.location)[-1] == self.cleavage_site:

                    try:
                        # qualifiers is structured like this:
                        # qualifiers: Key: note, Value: Cleavage; by factor Xa
                        note = feat.qualifiers["note"]
                        vals = note.split(";")

                        if vals[0] == "Cleavage":
                            protease = vals[1][4:]
                            proteases.append(protease)

                    except:
                        continue

        if len(proteases) >= 1:
            self.annot["protease_uniprot"] = "|".join(proteases)

    def merops_protease(self, merops, merops_name):
        """Queries MEROPS database static files for known annotated proteases
        cleaving in the cleavage site"""

        merops_acc = merops[merops.uniprot_acc == self.acc]
        proteases = []
        names = []

        for row in range(len(merops_acc)):

            # this is correct equality as p1 is amino acid before cleavage in actual
            # numbering, while cleavage_site is p1' index in python 0-index numbering
            if merops_acc.iloc[row]["p1"] == self.cleavage_site:

                proteases.append(merops_acc.iloc[row]["code"])

        for p in proteases:

            prot_name_df = merops_name[merops_name.code == p]
            if len(prot_name_df) == 0:
                names.append("Unknown")
            else:
                names.append(prot_name_df.iloc[0]["name"])

        if len(proteases) >= 1:
            self.annot["protease_merops_code"] = "|".join(proteases)
            self.annot["protease_merops_name"] = "|".join(names)
