import datetime
import json

from ord_schema.message_helpers import write_message
from ord_schema.proto.dataset_pb2 import *
from ord_schema.proto.reaction_pb2 import *

DATASET_DESCRIPTION = """C-H activation of oxazoles from the second demo case in the PyParser paper. \
It is also referred to as the "reaction optimization" campaign. \
Duplicate conditions are found in wells A2/A8, A3/A9, A5/A11 and A6/A12, then again in wells \
C2/C8, C3/C9, C5/C11 and C6/C12. This duplication was done to ensure reproducibility between the \
two 24-well metal heating blocks."""

SMILES_INTERNAL_STD = "c1ccc(CN(Cc2ccccc2)c2ccccc2)cc1"
SMILES_2 = "Brc1cnc2ccccc2c1"  # this is also the limiting reactant
SMILES_3 = "OCc1cocn1"
SMILES_4DOT1 = "OCc1coc(-c2cnc3ccccc3c2)n1"
SMILES_4DOT2 = "OCc1ncoc1-c1cnc2ccccc2c1"
SMILES_4DOT3 = "OCc1nc(-c2cnc3ccccc3c2)oc1-c1cnc2ccccc2c1"
SMILES_pd_catalyst_solvent = "O1C(C)CCC1"
SMILES_cesium_carbonate = "C(=O)([O-])[O-].[Cs+].[Cs+]"
SMILES_potassium_pivalate = "CC(C)(C)C(=O)[O-].[K+]"
SMILES_solvent = "C1COCCO1"
SMILES_copper_cocatalyst = "[Cu]Br.CSC"

with open("ocr/catalysts.json", "r") as f:
    CATALYST_INFO = json.load(f)


def build_reaction_smiles(catalyst_smiles: str, base_smiles: str) -> str:
    return (f"{SMILES_2}.{SMILES_3}>"
            f"{catalyst_smiles}.{base_smiles}.{SMILES_solvent}.{SMILES_copper_cocatalyst}"
            f">{SMILES_4DOT1}.{SMILES_4DOT2}.{SMILES_4DOT3}")


def build_compound(smi: str, amount: Amount, role: ReactionRole.ReactionRoleType, name: str = None, is_limiting=False):
    # TODO `ord_schema.message_helpers.build_compound` does not allow me to use an `Amount` instance
    cis = [CompoundIdentifier(
        type=CompoundIdentifier.CompoundIdentifierType.SMILES,
        value=smi,
    ), ]
    if name:
        ci = CompoundIdentifier(
            type=CompoundIdentifier.CompoundIdentifierType.NAME,
            value=name
        )
        cis.append(ci)
    return Compound(identifiers=cis, amount=amount, reaction_role=role, is_limiting=is_limiting)


def build_reaction(plate_row: str, plate_col_idx: int, corrp_std_4dot1: float, corrp_std_4dot2: float) -> Reaction:
    plate_row_idx = ["", "A", "B", "C", "D"].index(plate_row)  # 1-indexed

    if plate_row_idx <= 2:
        base_smiles = SMILES_cesium_carbonate
    else:
        base_smiles = SMILES_potassium_pivalate

    catalyst = CATALYST_INFO[((plate_row_idx - 1) * 12 + plate_col_idx - 1) % 24]
    catalyst_name = catalyst['name']
    catalyst_smi = catalyst['smiles']
    catalyst_cas = catalyst['cas']
    catalyst_mdl = catalyst['mdl']

    # create identifiers - reaction smiles, file names
    reaction_identifiers = [
        ReactionIdentifier(
            type=ReactionIdentifier.ReactionIdentifierType.REACTION_SMILES,
            value=build_reaction_smiles(catalyst_smi, base_smiles),
            is_mapped=False,
        ),
        ReactionIdentifier(
            type=ReactionIdentifier.ReactionIdentifierType.CUSTOM,
            details="the vial label on the 48-vial plate",
            value=f"{plate_row}-{plate_col_idx:02}"
        ),
    ]

    # create inputs - SMILES, quantities
    """
    Inside a glovebox, a solution of each palladium precatalyst was prepared by dissolving 20 μmol of each
    in a sufficient quantity of 2-methyltetrahydrofuran to give 2000 μL of solution. 100 μL of each solution
    was then dosed to 8 x 24 well plates, and the solvent was removed in vacuo using the Genevac for 1 h.
    The plates were stored in the glovebox and used within 1 year of creation.
    """
    # TODO how do we represent this type of addition?
    reaction_input_pd_catalyst = ReactionInput(
        components=[
            build_compound(
                smi=catalyst_smi,
                name=catalyst_name,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MICROMOLE, value=1)),
                role=ReactionRole.ReactionRoleType.CATALYST,
            ),
        ],
        cat_smiles = ReactionInput.component[0].identifier[2](type= "CAS_NUMBER"),
        addition_order=1,
        texture=Texture(type=Texture.TextureType.SOLID),
        addition_device=ReactionInput.AdditionDevice(type= 'CUSTOM', details= 'Palladium precatalyst (20 umol was dissolved in 2-methyltetrahydrofuran and diluted to 2 mL. 100 uL of precatalyst solution was dispensed into the vial, and then the solvent was removed in vacuo using the Genevac for 1 h.')
    )


    """Inside a glovebox, a plate of pre-dispensed palladium precatalysts (1 μmol) was taken and tumble \
    stirrer bars were added to each vial. Cesium carbonate (6.5 mg, 20 μmol, 2 eq.) was added by solid \
    addition (aided by a metal transfer template) to rows A and B of the two pre-dispensed plates, and \
    potassium pivalate (2.81 mg, 20 μmol, 2 eq.) was added to rows C and D."""
    # TODO technically, stirrer was has addition_order=2...
    reaction_input_base = ReactionInput(
        components=[
            build_compound(
                smi=base_smiles,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MICROMOLE, value=20)),
                role=ReactionRole.ReactionRoleType.REAGENT,
            ),
        ],
        addition_order=2,
        texture=Texture(type=Texture.TextureType.SOLID),
    )

    reaction_input_2 = ReactionInput(
        components=[
            build_compound(
                smi=SMILES_2,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MICROMOLE, value=10)),
                role=ReactionRole.ReactionRoleType.REACTANT,
                is_limiting=True,
            ),
            build_compound(
                smi=SMILES_solvent,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=25),
                              volume_includes_solutes=True),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_order=3,
        texture=Texture(type=Texture.TextureType.LIQUID),
    )

    reaction_input_3 = ReactionInput(
        components=[
            build_compound(
                smi=SMILES_3,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MICROMOLE, value=11)),
                role=ReactionRole.ReactionRoleType.REACTANT,
            ),
            build_compound(
                smi=SMILES_solvent,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=25),
                              volume_includes_solutes=True),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_order=4,
        texture=Texture(type=Texture.TextureType.LIQUID),
    )

    # TODO cocatalyst?
    reaction_input_cocatalyst = ReactionInput(
        components=[
            build_compound(
                smi=SMILES_copper_cocatalyst,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MICROMOLE, value=1)),
                role=ReactionRole.ReactionRoleType.CATALYST,
            ),
            build_compound(
                smi=SMILES_solvent,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=50),
                              volume_includes_solutes=True),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_order=5,
        texture=Texture(type=Texture.TextureType.LIQUID),
    )

    # reaction setup
    reaction_setup = ReactionSetup(
        vessel=Vessel(
            type=Vessel.VesselType.VIAL,
        ),
        is_automated=True,
        environment=ReactionSetup.ReactionEnvironment(
            type=ReactionSetup.ReactionEnvironment.ReactionEnvironmentType.CUSTOM,
            details="""
            Sealed 24-well metal heating blocks
            """
        ),
    )
    # TODO how to represent preparation mix in glove box?
    # TODO how to represent sealed vial?

    # reaction conditions
    reaction_conditions = ReactionConditions(
        temperature=TemperatureConditions(
            setpoint=Temperature(
                value=100,
                units=Temperature.TemperatureUnit.CELSIUS
            ),
            control=TemperatureConditions.TemperatureControl(
                type=TemperatureConditions.TemperatureControl.TemperatureControlType.CUSTOM,
                details="heated tumble stirrer"
            )
        ),
        pressure=PressureConditions(
            control=PressureConditions.PressureControl(
                type=PressureConditions.PressureControl.PressureControlType.SEALED,
                details="Sealed 24-well metal heating blocks"
            )
        ),
        conditions_are_dynamic=False,
    )

    # workups
    reaction_workups = [
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.TEMPERATURE,
            details="""The reactions were then allowed to cool to room temperature""",
            temperature=TemperatureConditions(
                setpoint=Temperature(value=25, units=Temperature.TemperatureUnit.CELSIUS)),
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ADDITION,
            details="""each reaction vial was diluted with 500 μL of a pre-prepared solution 
            consisting of 1 μmol dibenzylaniline (labelled as an internal standard for the purposes of PyParse analysis)
             and 20 μmol dimethylsulfoxide in acetonitrile.""",
            input=ReactionInput(
                components=[
                    build_compound(
                        smi=SMILES_INTERNAL_STD,
                        amount=Amount(moles=Moles(units=Moles.MolesUnit.MICROMOLE, value=1)),
                        role=ReactionRole.ReactionRoleType.INTERNAL_STANDARD,
                    ),
                    build_compound(
                        smi="CS(=O)C",
                        amount=Amount(moles=Moles(units=Moles.MolesUnit.MICROMOLE, value=20)),
                        role=ReactionRole.ReactionRoleType.WORKUP,
                    ),
                    build_compound(
                        smi="CCN",
                        amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=500),
                                      volume_includes_solutes=True),
                        role=ReactionRole.ReactionRoleType.WORKUP,
                    ),
                ],
                texture=Texture(type=Texture.TextureType.LIQUID),
            ),
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.STIRRING,
            details="""The reactions were sealed and shaken to ensure a homogenous mixture""",
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ALIQUOT,
            details="""50 μL of each reaction mixture was then aspirated and dispensed to the analysis plate""",
            amount=Amount(volume=Volume(value=50, units=Volume.VolumeUnit.MICROLITER))
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ADDITION,
            details="""Each well of the analysis plate was diluted with 150 μL of dimethylsulfoxide, and the analysis 
            plate was sealed with polypropylene film""",
            input=ReactionInput(
                components=[
                    build_compound(
                        smi="CS(=O)C",
                        amount=Amount(volume=Volume(value=150, units=Volume.VolumeUnit.MICROLITER)),
                        role=ReactionRole.ReactionRoleType.WORKUP,
                    ),
                ],
                texture=Texture(type=Texture.TextureType.LIQUID)
            )
        ),
        # TODO workup for "seal"?
    ]

    # outcomes
    # TODO relative peak area against internal standard
    # TODO there are typos in the subsection titles above figure S7 and figure S8
    outcome = ReactionOutcome(
        reaction_time=Time(value=18, units=Time.TimeUnit.HOUR),
        products=[
            ProductCompound(
                identifiers=[
                    CompoundIdentifier(type="SMILES", value=SMILES_4DOT1),
                ],
                is_desired_product=True,
                measurements=[
                    ProductMeasurement(
                        analysis_key="LC-MS",
                        type=ProductMeasurement.ProductMeasurementType.AREA,
                        details="""the metric used for analysis was “corrP/STD”; this metric was selected as it shades 
                        each well in the heatmap visualisation according to the normalised ratio of product to 
                        internal standard.""",
                        uses_internal_standard=True,
                        retention_time=Time(value=0.71, units=Time.TimeUnit.MINUTE),
                        percentage=Percentage(value=corrp_std_4dot1),
                    ),
                ],
                reaction_role=ReactionRole.ReactionRoleType.PRODUCT
            ),
            ProductCompound(
                identifiers=[
                    CompoundIdentifier(type="SMILES", value=SMILES_4DOT2),
                ],
                is_desired_product=False,
                measurements=[
                    ProductMeasurement(
                        analysis_key="LC-MS",
                        type=ProductMeasurement.ProductMeasurementType.AREA,
                        details="""the metric used for analysis was “corrP/STD”; this metric was selected as it shades 
                each well in the heatmap visualisation according to the normalised ratio of product to 
                internal standard.""",
                        uses_internal_standard=True,
                        retention_time=Time(value=0.73, units=Time.TimeUnit.MINUTE),
                        percentage=Percentage(value=corrp_std_4dot2),
                    ),
                ],
                reaction_role=ReactionRole.ReactionRoleType.SIDE_PRODUCT
            ),
            ProductCompound(
                identifiers=[
                    CompoundIdentifier(type="SMILES", value=SMILES_4DOT3),
                ],
                is_desired_product=False,
                measurements=[
                    ProductMeasurement(
                        analysis_key="LC-MS",
                        type=ProductMeasurement.ProductMeasurementType.CUSTOM,
                        details="only retention time is reported",
                        retention_time=Time(value=1.01, units=Time.TimeUnit.MINUTE),
                    ),
                ],
                reaction_role=ReactionRole.ReactionRoleType.SIDE_PRODUCT
            ),

        ],
        analyses={
            "LC-MS": Analysis(
                type=Analysis.AnalysisType.LCMS,
                details="""
                LCMS System B
                Column: 50 mm x 2.1 mm ID, Acquity UPLC CSH C18 column.
                Flow Rate: 1 mL/min
                Temp: 40 ˚C
                UV detection range: 210 to 350 nm
                Mass spectrum: Recorded on a mass spectrometer using alternate-scan positive and negative
                mode electrospray ionisation, with a mass range of 100 – 1000.
                Solvents:
                A: 10 mM solution of ammonium bicarbonate in water adjusted to pH 10 with ammonia solution.
                B: Acetonitrile
                Gradient: Time (min.) A% B%
                0 100 0
                1.5 3 97
                1.9 3 97
                2.0 100 0
                """,
            )
        }
    )

    notes = ReactionNotes(
        procedure_details="""
        - reaction input `pd_catalyst` was added as a solution (2-methyltetrahydrofuran 1 mmol/L) 
          and then the solvent was removed to leave solid solute before adding other reactants/reagents.
        - stirrer was added after `pd_catalyst` and before `base`
        - reaction inputs were prepared inside a glove box but the reaction itself did not happen in a glove box
        - the reaction vial was sealed after adding all reaction inputs
        """
    )

    reaction = Reaction(
        identifiers=reaction_identifiers,
        inputs={
            "2": reaction_input_2,
            "3": reaction_input_3,
            "pd_catalyst": reaction_input_pd_catalyst,
            "cocatalyst": reaction_input_cocatalyst,
            "base": reaction_input_base,
        },
        setup=reaction_setup,
        notes=notes,
        conditions=reaction_conditions,
        workups=reaction_workups,
        outcomes=[outcome],
        provenance=ReactionProvenance(
            doi="10.1039/d3dd00167a",
            record_created=RecordEvent(
                time=DateTime(value=datetime.datetime.now(tz=datetime.UTC).strftime("%Y-%m-%d %H:%M:%S %Z")),
                person=Person(username="qai", name="Qianxiang Ai", email="qai@mit.edu")),
        )
    )
    return reaction

if __name__ == '__main__':
    # figure S7, retention time 0.71
    CORRP_STD_4DOT1 = dict()
    CORRP_STD_4DOT1["A"] = [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0]
    CORRP_STD_4DOT1["B"] = [0, 2, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0]
    CORRP_STD_4DOT1["C"] = [31, 73, 6, 27, 38, 4, 10, 34, 5, 3, 24, 1]
    CORRP_STD_4DOT1["D"] = [16, 63, 80, 2, 48, 100, 56, 53, 24, 1, 13, 55]

    # figure S8, retention time 0.73
    CORRP_STD_4DOT2 = dict()
    CORRP_STD_4DOT2["A"] = [0, 0, 4, 1, 5, 4, 2, 0, 0, 6, 9, 1]
    CORRP_STD_4DOT2["B"] = [1, 11, 11, 32, 4, 7, 0, 0, 1, 20, 5, 7]
    CORRP_STD_4DOT2["C"] = [6, 0, 46, 8, 7, 74, 21, 6, 21, 17, 11, 100]
    CORRP_STD_4DOT2["D"] = [33, 2, 0, 79, 3, 0, 2, 0, 11, 6, 20, 0]

    REACTIONS = []
    for PLATE_ROW in ["A", "B", "C", "D"]:
        for PLATE_COL_ID in range(1, 13):
            REACTION = build_reaction(PLATE_ROW, PLATE_COL_ID, CORRP_STD_4DOT1[PLATE_ROW][PLATE_COL_ID - 1],
                                      CORRP_STD_4DOT2[PLATE_ROW][PLATE_COL_ID - 1])
            REACTIONS.append(REACTION)
    DATASET = Dataset(
        description=DATASET_DESCRIPTION,
        name="C-H activation of oxazoles from the second demo case in the PyParser paper",
        reactions=REACTIONS,
    )
    write_message(DATASET, "dataset.pbtxt")
