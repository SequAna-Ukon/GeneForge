#!/usr/bin/env python

# Memory-safe replacement for funannotate's aux_scripts/iprscan2annotations.py
# Original script written for funannotate by Jon Palmer (2017).
#
# BUG IN ORIGINAL: uses etree.iterparse() but never calls elem.clear(), so the
# entire XML tree is built and held in memory for the full duration of parsing.
# On large genomes (tens of thousands of proteins, multi-GB InterProScan XML),
# this grows unbounded and gets OOM-killed by the kernel with no traceback.
#
# FIX: single-pass streaming parse. Namespace stripping and data extraction
# happen in the same pass, per-<protein> element, and each <protein> element
# is cleared immediately after it's written to output. Memory stays roughly
# constant regardless of input file size.
#
# Usage identical to the original: iprscan2annotations.py IPRSCAN.xml OUTPUT.txt
# Requires FUNANNOTATE_DB env var to be set (same requirement as original).

import gzip
import sys
import os
import xml.etree.cElementTree as etree
from goatools import obo_parser


def convertGOattribute(namespacein):
    namespace = namespacein.upper()
    if namespace == "BIOLOGICAL_PROCESS":
        attribute = "go_process"
    elif namespace == "MOLECULAR_FUNCTION":
        attribute = "go_function"
    elif namespace == "CELLULAR_COMPONENT":
        attribute = "go_component"
    else:
        attribute = "go_unknown"
    return attribute


def strip_ns(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


def main():
    """Streaming, memory-safe version of interpro annotations to tab delimited script."""

    if len(sys.argv) < 3:
        print("Usage: iprscan2annotations_streaming.py IPRSCAN.xml OUTPUT.annotations.txt")
        sys.exit(1)

    goDict = {}
    for item in obo_parser.OBOReader(
        os.path.join(os.environ["FUNANNOTATE_DB"], "go.obo")
    ):
        namespace = convertGOattribute(item.namespace)
        goDict[item.id] = {"name": item.name, "namespace": namespace}
        for nm in item.alt_ids:
            goDict[nm] = {"name": item.name, "namespace": namespace}

    _opener = gzip.open(sys.argv[1], "rt") if sys.argv[1].endswith(".gz") else open(sys.argv[1])
    with open(sys.argv[2], "w") as output:
        with _opener as xml_file:
            context = etree.iterparse(xml_file, events=("end",))
            for _, elem in context:
                # strip namespace from this element's tag and attributes,
                # same as the original script's first pass, but done inline
                # per-element as we stream instead of over the whole tree.
                elem.tag = strip_ns(elem.tag)
                for at in list(elem.attrib.keys()):
                    if "}" in at:
                        newat = strip_ns(at)
                        elem.attrib[newat] = elem.attrib[at]
                        del elem.attrib[at]

                if elem.tag != "protein":
                    continue

                hits = elem
                IDs = []
                iprs = []
                gos = {}
                signalp = []
                for lv1 in hits:
                    if lv1.tag == "xref":
                        name = lv1.get("id")
                        IDs.append(name)
                    if lv1.tag == "matches":
                        for e in lv1.findall(".//entry"):
                            if e.get("ac") not in iprs:
                                iprs.append(e.get("ac"))
                        for g in lv1.findall(".//go-xref"):
                            cat = g.get("category", None)
                            goID = g.get("id", None)
                            desc = g.get("name", None)
                            if not goID:
                                continue
                            if not cat or not desc:
                                if goID in goDict:
                                    cat = goDict[goID]["namespace"]
                                    desc = goDict[goID]["name"]
                                else:
                                    continue
                            else:
                                cat = convertGOattribute(cat)
                            goHit = (cat, desc, goID)
                            if goID not in gos:
                                gos[goID] = goHit
                        for s in lv1.findall(".//signalp-match"):
                            for lib in s.findall(".//signature-library-release"):
                                if lib.get("library") == "SIGNALP_EUK":
                                    for loc in s.findall(".//signalp-location"):
                                        signalp.append(
                                            (loc.get("start"), loc.get("end"))
                                        )

                if len(iprs) > 0:
                    for i in IDs:
                        for x in iprs:
                            output.write(f"{i}\tdb_xref\tInterPro:{x}\n")
                if len(gos) > 0:
                    for i in IDs:
                        for goid in gos:
                            x = gos[goid]
                            GOID = x[2].replace("GO:", "")
                            output.write(f"{i}\t{x[0]}\t{x[1]}|{GOID}||IEA\n")

                # release this protein's subtree now that it's written out
                elem.clear()


if __name__ == "__main__":
    main()

