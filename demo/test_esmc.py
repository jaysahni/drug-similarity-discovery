"""Does ESM-C actually separate related proteins from unrelated ones?

A test that only checks `embed()` returns the right shape would pass on random
numbers. So this checks the embeddings carry biology, on sequences whose
relationships are known independently of the model:

  * four human CDKs (CDK1/2/4/6) -- one family, one fold, homologous
  * four unrelated human proteins -- insulin, lysozyme, albumin, haemoglobin beta
  * CDK2 with its residues shuffled -- same amino acid composition, no protein

The orderings asserted below were written before the numbers were seen, and none
of them is a tuned threshold. If ESM-C fails to separate these, the test fails
and that is the finding; it does not get relaxed until it passes.

Every sequence is the full UniProt canonical entry, accession given inline, so
the inputs are checkable rather than paraphrased.

Run:  ./env-kit/bin/python -m demo.test_esmc
"""

from __future__ import annotations

import random
import time

import numpy as np
from scipy.stats import mannwhitneyu

from demo import esmc

SEQ: dict[str, str] = {
    # CDK2_HUMAN, UniProt P24941, 298 aa
    "CDK2": (
        "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKL"
        "LDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQN"
        "LLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVT"
        "RRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLS"
        "QMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL"
    ),
    # CDK6_HUMAN, UniProt Q00534, 326 aa
    "CDK6": (
        "MEKDGLCRADQQYECVAEIGEGAYGKVFKARDLKNGGRFVALKRVRVQTGEEGMPLSTIREVAVLR"
        "HLETFEHPNVVRLFDVCTVSRTDRETKLTLVFEHVDQDLTTYLDKVPEPGVPTETIKDMMFQLLRG"
        "LDFLHSHRVVHRDLKPQNILVTSSGQIKLADFGLARIYSFQMALTSVVVTLWYRAPEVLLQSSYAT"
        "PVDLWSVGCIFAEMFRRKPLFRGSSDVDQLGKILDVIGLPGEEDWPRDVALPRQAFHSKSAQPIEK"
        "FVTDIDELGKDLLLKCLTFNPAKRISAYSALSHPYFQDLERCKENLDSHLPPSQNTSELNTA"
    ),
    # CDK1_HUMAN, UniProt P06493, 297 aa
    "CDK1": (
        "MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRHPNIVSL"
        "QDVLMQDSRLYLIFEFLSMDLKKYLDSIPPGQYMDSSLVKSYLYQILQGIVFCHSRRVLHRDLKPQ"
        "NLLIDDKGTIKLADFGLARAFGIPIRVYTHEVVTLWYRSPEVLLGSARYSTPVDIWSIGTIFAELA"
        "TKKPLFHGDSEIDQLFRIFRALGTPNNEVWPEVESLQDYKNTFPKWKPGSLASHVKNLDENGLDLL"
        "SKMLIYDPAKRISGKMALNHPYFNDLDNQIKKM"
    ),
    # CDK4_HUMAN, UniProt P11802, 303 aa
    "CDK4": (
        "MATSRYEPVAEIGVGAYGTVYKARDPHSGHFVALKSVRVPNGGGGGGGLPISTVREVALLRRLEAF"
        "EHPNVVRLMDVCATSRTDREIKVTLVFEHVDQDLRTYLDKAPPPGLPAETIKDLMRQFLRGLDFLH"
        "ANCIVHRDLKPENILVTSGGTVKLADFGLARIYSYQMALTPVVVTLWYRAPEVLLQSTYATPVDMW"
        "SVGCIFAEMFRRKPLFCGNSEADQLGKIFDLIGLPPEDDWPRDVSLPRGAFPPRGPRPVQSVVPEM"
        "EESGAQLLLEMLTFNPHKRISAFRALQHSYLHKDEGNPE"
    ),
    # INS_HUMAN, UniProt P01308, 110 aa
    "INS": (
        "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQV"
        "ELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    ),
    # LYZ_HUMAN, UniProt P61626, 148 aa
    "LYZ": (
        "MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAG"
        "DRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRN"
        "RCQNRDVRQYVQGCGV"
    ),
    # ALB_HUMAN, UniProt P02768, 609 aa
    "ALB": (
        "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKL"
        "VNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDD"
        "NPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADK"
        "AACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTK"
        "VHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLA"
        "ADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAK"
        "VFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCC"
        "KHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNA"
        "ETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAE"
        "EGKKLVAASQAALGL"
    ),
    # HBB_HUMAN, UniProt P68871, 147 aa
    "HBB": (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGK"
        "KVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQ"
        "KVVAGVANALAHKYH"
    ),
    # HBA_HUMAN, UniProt P69905, 142 aa
    "HBA": (
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADA"
        "LTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLAS"
        "VSTVLTSKYR"
    ),
    # IGHG1_HUMAN, UniProt P01857, 399 aa
    "IGHG1": (
        "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLS"
        "SVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKD"
        "TLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLN"
        "GKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEW"
        "ESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPEL"
        "QLEESCAEAQDGELDGLWTTITIFITLFLLSVCYSATVTFFKVKWIFSSVVDLKQTIIPDYRNMIG"
        "QGA"
    ),
    # IGKC_HUMAN, UniProt P01834, 107 aa
    "IGKC": (
        "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTY"
        "SLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
    ),}

KINASES = ["CDK1", "CDK2", "CDK4", "CDK6"]
OTHERS = ["INS", "LYZ", "ALB", "HBB"]

failures: list[str] = []


def check(name: str, ok: bool, detail: str) -> None:
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: {detail}")
    if not ok:
        failures.append(f"{name}: {detail}")


def test_identity(vectors: dict[str, np.ndarray]) -> None:
    """A sequence against itself must be cosine 1."""
    print("\n1. identity")
    worst = max(
        (abs(esmc.similarity(v, v) - 1.0), k) for k, v in vectors.items()
    )
    check(
        "self-cosine == 1",
        worst[0] < 1e-6,
        f"largest deviation {worst[0]:.2e} (on {worst[1]}), n={len(vectors)}",
    )

    # Re-embedding the same string must hit the cache and return the identical
    # vector; a cache keyed wrongly would show up here as a nonzero difference.
    again = esmc.embed([SEQ["CDK2"]])[0]
    check(
        "cache round-trip is exact",
        float(np.abs(again - vectors["CDK2"]).max()) == 0.0,
        f"max |diff| {float(np.abs(again - vectors['CDK2']).max()):.2e}",
    )


def test_family(vectors: dict[str, np.ndarray]) -> None:
    """Within-family cosines should sit above between-family cosines."""
    print("\n2. CDK family vs unrelated proteins")
    within, between = [], []
    for i, a in enumerate(KINASES):
        for b in KINASES[i + 1 :]:
            within.append((f"{a}~{b}", esmc.similarity(vectors[a], vectors[b])))
        for b in OTHERS:
            between.append((f"{a}~{b}", esmc.similarity(vectors[a], vectors[b])))

    for label, value in within:
        print(f"      within  {label:<12} {value:.4f}")
    for label, value in between:
        print(f"      between {label:<12} {value:.4f}")

    w = np.array([v for _, v in within])
    b = np.array([v for _, v in between])
    u, p = mannwhitneyu(w, b, alternative="greater")

    check(
        "every within-family pair beats every between-family pair",
        w.min() > b.max(),
        f"min within {w.min():.4f} ({min(within, key=lambda t: t[1])[0]}) vs "
        f"max between {b.max():.4f} ({max(between, key=lambda t: t[1])[0]})",
    )
    check(
        "within > between (Mann-Whitney U, one-sided)",
        p < 0.05,
        f"within mean {w.mean():.4f} (n={len(w)}), between mean {b.mean():.4f} "
        f"(n={len(b)}), U={u:.0f}, p={p:.2e}",
    )


def test_shuffle(vectors: dict[str, np.ndarray]) -> None:
    """A composition-matched shuffle is not a protein, and should score below a homolog."""
    print("\n3. shuffled control")
    rng = random.Random(0)
    residues = list(SEQ["CDK2"])
    rng.shuffle(residues)
    shuffled = "".join(residues)
    assert sorted(shuffled) == sorted(SEQ["CDK2"]), "shuffle changed composition"

    v_shuf = esmc.embed([shuffled])[0]
    to_shuffle = esmc.similarity(vectors["CDK2"], v_shuf)
    to_homolog = esmc.similarity(vectors["CDK2"], vectors["CDK6"])
    to_stranger = esmc.similarity(vectors["CDK2"], vectors["INS"])

    print(f"      CDK2 ~ shuffled CDK2  {to_shuffle:.4f}  (same composition, seed 0)")
    print(f"      CDK2 ~ CDK6           {to_homolog:.4f}  (homolog)")
    print(f"      CDK2 ~ insulin        {to_stranger:.4f}  (unrelated)")
    check(
        "homolog beats shuffle",
        to_homolog > to_shuffle,
        f"{to_homolog:.4f} > {to_shuffle:.4f}, margin {to_homolog - to_shuffle:+.4f}",
    )


def test_multichain(vectors: dict[str, np.ndarray]) -> None:
    """Two-chain entities: all four combine modes, on real two-chain assemblies."""
    print("\n4. multi-chain entities")
    entities = {
        "haemoglobin (HBA+HBB)": ["HBA", "HBB"],
        "IgG1 constant (IGHG1+IGKC)": ["IGHG1", "IGKC"],
    }
    for label, parts in entities.items():
        print(f"      {label}")
        for mode in esmc.COMBINE:
            v = esmc.embed_entity([SEQ[p] for p in parts], combine=mode)
            if mode == "concat":
                print(f"        {mode:<11} dim {v.shape[0]} (not in the single-chain space)")
                continue
            own = [esmc.similarity(v, vectors[p]) for p in parts]
            far = esmc.similarity(v, vectors["ALB"])
            print(
                f"        {mode:<11} dim {v.shape[0]}  "
                + "  ".join(f"~{p} {s:.4f}" for p, s in zip(parts, own))
                + f"  ~ALB {far:.4f}"
            )
            check(
                f"{label} [{mode}] closer to its own chains than to albumin",
                min(own) > far,
                f"min own {min(own):.4f} > {far:.4f}",
            )

    # combine='concat' must refuse to stack entities of differing chain count.
    try:
        esmc.embed_entities([[SEQ["HBA"], SEQ["HBB"]], [SEQ["ALB"]]], combine="concat")
        check("concat refuses ragged chain counts", False, "it did not raise")
    except ValueError as exc:
        check("concat refuses ragged chain counts", True, str(exc).split(" --")[0])


def test_limits() -> None:
    """Long sequences and bad alphabets must fail loudly."""
    print("\n5. guards")
    over = (SEQ["ALB"] * 10)[: esmc.MAX_RESIDUES + 1]
    try:
        esmc.embed([over])
        check("over-context sequence refused", False, "it did not raise")
    except ValueError as exc:
        check(
            "over-context sequence refused",
            True,
            f"{len(over)} aa > {esmc.MAX_RESIDUES}: {str(exc).split('.')[0]}",
        )
    check(
        "longest real chain here fits",
        max(len(s) for s in SEQ.values()) <= esmc.MAX_RESIDUES,
        f"longest is ALB at {len(SEQ['ALB'])} aa, limit {esmc.MAX_RESIDUES}",
    )
    try:
        esmc.embed([SEQ["INS"].lower()])
        check("lowercase sequence refused", False, "it did not raise")
    except ValueError:
        check("lowercase sequence refused", True, "raised ValueError as intended")


def main() -> None:
    names = KINASES + OTHERS + ["HBA", "IGHG1", "IGKC"]
    print(f"model {esmc.MODEL} on {esmc.device()}")

    t0 = time.perf_counter()
    esmc.load()
    load_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    matrix = esmc.embed([SEQ[n] for n in names])
    embed_s = time.perf_counter() - t0

    residues = sum(len(SEQ[n]) for n in names)
    print(
        f"loaded in {load_s:.2f}s; embedded n={len(names)} sequences "
        f"({residues} residues) in {embed_s:.2f}s "
        f"= {embed_s / len(names) * 1000:.0f} ms/sequence, dim {matrix.shape[1]}"
    )
    vectors = dict(zip(names, matrix))

    test_identity(vectors)
    test_family(vectors)
    test_shuffle(vectors)
    test_multichain(vectors)
    test_limits()

    print()
    if failures:
        print(f"{len(failures)} CHECK(S) FAILED -- this is a result, not a bug to hide:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
