import datetime
import json

import pandas as pd
from ord_schema.message_helpers import write_message
from ord_schema.proto.dataset_pb2 import *
from ord_schema.proto.reaction_pb2 import *
from rdkit.Chem import AllChem, MolFromSmiles, MolToSmiles

with open("mappings.json", "r") as f:
    MAPPINGS = json.load(f)

MAPPING_P2NAME, MAPPING_P2SMILES, MAPPING_L2SMILES, MAPPING_S2SMILES = MAPPINGS
COMPOUND_3B_SMILES = "C#CCCCCC=O"
COMPOUND_BASE_SMILES = "Cc1cccc(C)n1"  # 2,6-Lutidine
COMPOUND_SOLVENT_SMILES = "CN(C)C=O"  # DMF
COMPOUND_INTERNAL_STANDARD_SMILES = "OC(C#C)c1ccccc1"
COMPOUND_INTERNAL_STANDARD_NAME = "1-phenylprop-2-yn-1-ol"
COMPOUND_D3_SMILES = "[N-]=[N+]=Nc1ccc(C(=O)OC[C@H](Cc2ccccc2)NC(=O)OCC2c3ccccc3-c3ccccc32)cc1"  # derivatization reactant

REACTION_RECORDS = pd.read_csv('sheet0.csv').to_dict(orient='records')
REACTION_SMARTS_LHS = "[C:1]#[C:2][C:3][C:4][C:5][C:6][C:7]=[O:8].[c:9]1[c:10][c:11][c:12]([C:13](=[O:14])[C:15]Br)[c:17][c:18]1"
REACTION_SMARTS_RHS_S = "[c:9]1[c:10][c:11][c:12]([C:13](=[O:14])[C:15][C@@H:6]([C:5][C:4][C:3][C:2]#[C:1])([C:7]=[O:8]))[c:17][c:18]1"
REACTION_SMARTS_RHS_R = "[c:9]1[c:10][c:11][c:12]([C:13](=[O:14])[C:15][C@H:6]([C:5][C:4][C:3][C:2]#[C:1])([C:7]=[O:8]))[c:17][c:18]1"

REACTION_SMARTS_R = REACTION_SMARTS_LHS + ">>" + REACTION_SMARTS_RHS_R
REACTION_SMARTS_S = REACTION_SMARTS_LHS + ">>" + REACTION_SMARTS_RHS_S
RXN_R = AllChem.ReactionFromSmarts(REACTION_SMARTS_R)
RXN_S = AllChem.ReactionFromSmarts(REACTION_SMARTS_S)

REACTION_SMARTS_DERIVATIZATION = "[N-:1]=[N+:2]=[N:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[N:1]=[N:2]-[N:3]1"
RXN_DERIVATIZATION = AllChem.ReactionFromSmarts(REACTION_SMARTS_DERIVATIZATION)


def str2percentage(s: str) -> Percentage:
    assert s.endswith("%")
    if s.startswith("--"):
        s = s.replace("--", "")
    return Percentage(value=float(s.strip("%")))


def get_reaction_smiles(p_smi, s_smi, l_smi):
    product_s_mol = RXN_S.RunReactants([MolFromSmiles(COMPOUND_3B_SMILES), MolFromSmiles(s_smi)])[0][0]
    product_s_smi = MolToSmiles(product_s_mol)
    product_r_mol = RXN_R.RunReactants([MolFromSmiles(COMPOUND_3B_SMILES), MolFromSmiles(s_smi)])[0][0]
    product_r_smi = MolToSmiles(product_r_mol)
    return f"{COMPOUND_3B_SMILES}.{s_smi}>{l_smi}.{p_smi}.{COMPOUND_SOLVENT_SMILES}.{COMPOUND_BASE_SMILES}>{product_s_smi}.{product_r_smi}", product_s_smi, product_r_smi


def build_compound(smi: str, amount: Amount, role: ReactionRole.ReactionRoleType, name: str = None):
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
    return Compound(identifiers=cis, amount=amount, reaction_role=role, )


def create_reaction(reaction_record) -> Reaction:
    # parse reaction record
    S_label = reaction_record["Substrate"]  # Sx
    S_smiles = MAPPING_S2SMILES[S_label]

    P_label = reaction_record["Photocatalyst"]  # Px
    P_smiles = MAPPING_P2SMILES[P_label]
    P_name = MAPPING_P2NAME[P_label]

    L_label = reaction_record["Organocatalyst"]  # Lx
    L_smiles = MAPPING_L2SMILES[L_label]

    A_R = reaction_record["AR"]
    """
    > By adding another alkyne-containing compound 1-phenylprop-2-
    > yn-1-ol (4) as an internal standard, the relative MS yields of different
    > reaction systems could be compared by calculating the relative peak
    > areas of ions (AR = Aproduct/Ainternal standard)33 between the derivatives
    > from the reaction products and internal standard while determining the
    > ee values. 
    """

    ee = reaction_record["ee (IM-MS)"]
    """
    The reported *ee* value was calculated s.t. positive - R config, negative - S config
    > Among these catalysts, the chiral products obtained with L2 and L3 are mainly R-configuration, 
    > while those obtained with L4 and L10 are mainly S-configuration.
    """

    file_name = reaction_record["File name"]

    reaction_smiles, product_s_smi, product_r_smi = get_reaction_smiles(P_smiles, S_smiles, L_smiles)

    # create identifiers - reaction smiles, file names
    reaction_identifiers = [
        ReactionIdentifier(
            type=ReactionIdentifier.ReactionIdentifierType.REACTION_SMILES,
            details="two products are enantiomers",
            value=reaction_smiles,
            is_mapped=False,
            # TODO use mappers
        ),
        ReactionIdentifier(
            type=ReactionIdentifier.ReactionIdentifierType.CUSTOM,
            details="Each reaction in this dataset is assigned a unique `File name` partially describing the reactants",
            value=file_name
        ),
    ]

    # create inputs - SMILES, quantities, addition device
    addition_device = ReactionInput.AdditionDevice(
        type=ReactionInput.AdditionDevice.AdditionDeviceType.PIPETTE,
        details="VERSA 110, Aurora Biomed, Canada"
    )
    """
     we prepared the reaction components stock solutions as following: S1-S10 (2 M in DMF, 5 mL), 
     L1-L11 (400 mM in DMF, 5 mL), P1-P13 (10 mM in DMF, 5 mL), 3a and 2,6-lutidine mixture (4 M in DMF, 40 mL). 
     Then combinatorial mixed these components with equal 25 uL to 96-well plates by the automated liquid handing 
     system.
    """
    reaction_input_s = ReactionInput(
        components=[
            build_compound(
                smi=S_smiles,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MOLE, value=2 * 25 * 1e-6)),
                role=ReactionRole.ReactionRoleType.REACTANT,
            ),
            build_compound(
                smi=COMPOUND_SOLVENT_SMILES,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=25),
                              volume_includes_solutes=True),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_device=addition_device,
        addition_temperature=Temperature(value=298, units=Temperature.TemperatureUnit.KELVIN),
        # TODO how to describe room temperature?
        texture=Texture(type=Texture.TextureType.LIQUID)
    )

    reaction_input_l = ReactionInput(
        components=[
            build_compound(
                smi=L_smiles,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MOLE, value=0.4 * 25 * 1e-6)),
                role=ReactionRole.ReactionRoleType.CATALYST,
            ),
            build_compound(
                smi=COMPOUND_SOLVENT_SMILES,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=25)),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_device=addition_device,
        addition_temperature=Temperature(value=298, units=Temperature.TemperatureUnit.KELVIN),
        texture=Texture(type=Texture.TextureType.LIQUID)
    )

    reaction_input_p = ReactionInput(
        components=[
            build_compound(
                smi=P_smiles,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MOLE, value=10 * 1e-3 * 25 * 1e-6)),
                role=ReactionRole.ReactionRoleType.CATALYST,
                name=P_name,
            ),
            build_compound(
                smi=COMPOUND_SOLVENT_SMILES,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=25)),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_device=addition_device,
        addition_temperature=Temperature(value=298, units=Temperature.TemperatureUnit.KELVIN),
        texture=Texture(type=Texture.TextureType.LIQUID)
    )

    reaction_input_3b = ReactionInput(
        components=[
            build_compound(
                smi=COMPOUND_3B_SMILES,
                # TODO there seems to be a typo in methods, 6-Heptynal is 3a there while in Fig.3a it is 3b
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MOLE, value=4 * 25 * 1e-6)),
                # TODO it is a bit unclear if 4M is used for 3a or 2M is used
                role=ReactionRole.ReactionRoleType.REACTANT,
            ),
            build_compound(
                smi=COMPOUND_SOLVENT_SMILES,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=25)),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_device=addition_device,
        addition_temperature=Temperature(value=298, units=Temperature.TemperatureUnit.KELVIN),
        texture=Texture(type=Texture.TextureType.LIQUID)
    )

    reaction_input_base = ReactionInput(
        components=[
            build_compound(
                smi=COMPOUND_BASE_SMILES,
                amount=Amount(moles=Moles(units=Moles.MolesUnit.MOLE, value=4 * 25 * 1e-6)),
                role=ReactionRole.ReactionRoleType.REACTANT,
            ),
            build_compound(
                smi=COMPOUND_SOLVENT_SMILES,
                amount=Amount(volume=Volume(units=Volume.VolumeUnit.MICROLITER, value=25)),
                role=ReactionRole.ReactionRoleType.SOLVENT,
            ),
        ],
        addition_device=addition_device,
        addition_temperature=Temperature(value=298, units=Temperature.TemperatureUnit.KELVIN),
        texture=Texture(type=Texture.TextureType.LIQUID),
    )

    # reaction setup
    reaction_setup = ReactionSetup(
        vessel=Vessel(
            type=Vessel.VesselType.WELL_PLATE,
            details="450 uL Nunc U96 Microwell, Thermo Scientific",
            preparations=[
                VesselPreparation(
                    type=VesselPreparation.VesselPreparationType.PURGED,
                    details="A transparent acrylic top layer was fixed to the container with 12 flange bolts, "
                            "then the container was degassed with an oil pump for three times and refilled with "
                            "nitrogen and the container was connected with a nitrogen balloon to ensure a nitrogen "
                            "atmosphere for the reaction."
                ),
            ]
        ),
        is_automated=True,
    )

    # conditions
    reaction_conditions = ReactionConditions(
        temperature=TemperatureConditions(
            setpoint=Temperature(value=22, units=Temperature.TemperatureUnit.CELSIUS, precision=5),
            control=TemperatureConditions.TemperatureControl(
                type=TemperatureConditions.TemperatureControl.TemperatureControlType.CUSTOM,
                details=""" two 96-well plates were placed side by side into a home-made aluminum alloy container (see Supplementary Fig. 6). """
            )
        ),
        # TODO it is a bit conflicting that Fig. 3a is RT but there is an ice-cooler used in methods
        illumination=IlluminationConditions(
            type=IlluminationConditions.IlluminationType.LED,
            details="100 W LED compact fluorescent lamp (CFL) w",
            color="white",
            distance_to_vessel=Length(value=8, units=Length.LengthUnit.CENTIMETER)
        ),
        pressure=PressureConditions(
            control=PressureConditions.PressureControl(
                type=PressureConditions.PressureControl.PressureControlType.SEALED,
                details=""" each 96-well plate was covered with an optical glass to minimize the volatilization of solvent and components """
            ),
        ),
        stirring=StirringConditions(
            type=StirringConditions.StirringMethodType.SONICATION,
            details=""" After placing the 96-well plates in ultrasonic water bath for 10 seconds to mix the reaction uniformly, """,
        ),
        conditions_are_dynamic=False,
    )

    # workups -- the derivation reaction
    reaction_workups = [
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ALIQUOT,
            details="""4 uL mixture of each reacted solution was collected and diluted with 36 uL ACN""",
            amount=Amount(volume=Volume(value=4, units=Volume.VolumeUnit.MICROLITER))
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ADDITION,
            details="""4 uL mixture of each reacted solution was collected and diluted with 36 uL ACN""",
            input=ReactionInput(
                components=[
                    build_compound(
                        smi="CC#N",
                        amount=Amount(volume=Volume(value=36, units=Volume.VolumeUnit.MICROLITER)),
                        role=ReactionRole.WORKUP,
                    ),
                ],
                addition_device=addition_device,
            ),
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ADDITION,
            details="chiral resolving regent (150 mM in DCM, 40 uL)",
            input=ReactionInput(
                components=[
                    build_compound(
                        smi=COMPOUND_D3_SMILES,
                        amount=Amount(moles=Moles(value=150 * 1e-3 * 40 * 1e-6, units=Moles.MolesUnit.MOLE)),
                        role=ReactionRole.WORKUP,  # TODO if workup is a reaction what should the role be?
                    ),
                    build_compound(
                        smi=COMPOUND_SOLVENT_SMILES,
                        amount=Amount(volume=Volume(value=40, units=Volume.VolumeUnit.MICROLITER)),
                        role=ReactionRole.WORKUP,
                    ),
                ],
                addition_device=addition_device,
                # TODO ReactionWorkup should be sequential but ReactionInput here also has an addition_order
            ),
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ADDITION,
            details="CuI with DIPEA (400 mM in ACN) 40 uL",
            input=ReactionInput(
                components=[
                    build_compound(
                        smi="CCN(C(C)C)C(C)C.[Cu+].[I-]",
                        amount=Amount(moles=Moles(value=400 * 1e-3 * 40 * 1e-6, units=Moles.MolesUnit.MOLE)),
                        role=ReactionRole.WORKUP,  # TODO if workup is a reaction what should the role be?
                    ),
                    build_compound(
                        smi="CC#N",
                        amount=Amount(volume=Volume(value=40, units=Volume.VolumeUnit.MICROLITER)),
                        role=ReactionRole.WORKUP,
                    ),
                ],
                addition_device=addition_device,
            ),
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ADDITION,
            details="internal standard 1-phenylprop-2-yn-1-ol (4) (25 mM in ACN) 40 uL",
            input=ReactionInput(
                components=[
                    build_compound(
                        smi=COMPOUND_INTERNAL_STANDARD_SMILES,
                        amount=Amount(moles=Moles(value=25 * 1e-3 * 40 * 1e-6, units=Moles.MolesUnit.MOLE)),
                        role=ReactionRole.INTERNAL_STANDARD,
                        name=COMPOUND_INTERNAL_STANDARD_NAME
                    ),
                    build_compound(
                        smi="CC#N",
                        amount=Amount(volume=Volume(value=40, units=Volume.VolumeUnit.MICROLITER)),
                        role=ReactionRole.WORKUP,
                    ),
                ],
                addition_device=addition_device,
            ),
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.STIRRING,
            duration=Time(value=10, units=Time.TimeUnit.MINUTE),
            temperature=TemperatureConditions(
                setpoint=Temperature(value=298, units=Temperature.TemperatureUnit.KELVIN)),
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ALIQUOT,
            details="""After the reaction, 10 uL of supernatant for each reaction was transferred to another 96-well plate""",
            amount=Amount(volume=Volume(value=10, units=Volume.VolumeUnit.MICROLITER))
        ),
        ReactionWorkup(
            type=ReactionWorkup.ReactionWorkupType.ADDITION,
            details="After the reaction, 10 uL of supernatant for each reaction was transferred to another 96-well plate and diluted with 190 uL ACN before the IM-MS analysis",
            input=ReactionInput(
                components=[
                    build_compound(
                        smi="CC#N",
                        amount=Amount(volume=Volume(value=190, units=Volume.VolumeUnit.MICROLITER)),
                        role=ReactionRole.WORKUP,
                    ),
                ],
                addition_device=addition_device,
            ),
        ),
    ]

    # since peak area was reported for both enantiomers, we need to calculate individual enantiomer peak areas
    ee_r_percentage = str2percentage(ee)
    assert ee_r_percentage.value <= 100
    A_R_r = A_R * (1 + ee_r_percentage.value / 100) / 2
    A_R_s = A_R - A_R_r
    print(ee_r_percentage.value, A_R_r, A_R_s)

    # reaction outcomes
    product_compounds = []
    for product_smi, ee_string, peak_area, custom_identifier in zip(
            [product_r_smi, product_s_smi], [ee, "-" + ee], [A_R_r, A_R_s], ["R config", "S config"]
    ):
        product_compound = ProductCompound(
            identifiers=[
                CompoundIdentifier(type="SMILES", value=product_smi),
                CompoundIdentifier(type=CompoundIdentifier.CompoundIdentifierType.CUSTOM, value=custom_identifier,
                                   details="enantiomer label")
            ],
            is_desired_product=True,
            measurements=[
                ProductMeasurement(
                    analysis_key="IM-MS",
                    type=ProductMeasurement.ProductMeasurementType.SELECTIVITY,
                    uses_internal_standard=True,
                    percentage=str2percentage(ee_string),
                    selectivity=ProductMeasurement.Selectivity(
                        type=ProductMeasurement.Selectivity.SelectivityType.EE,
                        details=ee_string,
                    )
                ),
            ],
            reaction_role=ReactionRole.ReactionRoleType.PRODUCT
        )
        product_compounds.append(product_compound)

    for product_smi, ee_string, peak_area, custom_identifier in zip(
            [product_r_smi, product_s_smi], [ee, "-" + ee], [A_R_r, A_R_s],
            ["R config derivative", "S config derivative"]
    ):
        d_product = RXN_DERIVATIZATION.RunReactants((MolFromSmiles(COMPOUND_D3_SMILES), MolFromSmiles(product_smi)))[0][
            0]
        d_product_smi = MolToSmiles(d_product)
        # print(custom_identifier, d_product_smi)

        product_compound = ProductCompound(
            identifiers=[
                CompoundIdentifier(type="SMILES", value=d_product_smi),
                CompoundIdentifier(type=CompoundIdentifier.CompoundIdentifierType.CUSTOM, value=custom_identifier,
                                   details="enantiomer label")
            ],
            is_desired_product=False,
            measurements=[
                ProductMeasurement(
                    analysis_key="IM-MS",
                    type=ProductMeasurement.ProductMeasurementType.AREA,
                    uses_internal_standard=True,
                    details="""By adding another alkyne-containing compound 1-phenylprop-2-yn-1-ol (4) as an internal standard, the relative MS yields of different reaction systems could be compared by calculating the relative peakareas of ions (AR = Aproduct/Ainternal standard)""",
                    float_value=FloatValue(value=peak_area),
                )
            ],
            reaction_role=ReactionRole.ReactionRoleType.PRODUCT
        )
        product_compounds.append(product_compound)

    outcome = ReactionOutcome(
        reaction_time=Time(value=8, units=Time.TimeUnit.HOUR),
        products=product_compounds,
        analyses={
            "IM-MS": Analysis(
                type=Analysis.AnalysisType.MS,
                details="""Ion mobility (IM)-mass spectrometry (MS) analysis were performed on timsTOF Pro mass spectrometer (Bruker
Daltonics, Germany) with a CaptiveSpray ion source. The final solutions were directly injected into the ion source
by ultimate 3000 LC autosampler (Thermo Scientific, USA) with mobile phase. High-resolution (HR) MS data was
acquired with Orbitrap EliteTM mass spectrometer (Thermo Scientific, USA) with the following instrument
parameters: FTMS positive mode, spray voltage: 3.8 kV, source heater temperature: 350 oC, sheath gas flow rate: 40,
aux gas flow rate: 10."""
            )
        }
    )

    reaction = Reaction(
        identifiers=reaction_identifiers,
        inputs={
            "S": reaction_input_s,
            "L": reaction_input_l,
            "P": reaction_input_p,
            "3b": reaction_input_3b,
            "base": reaction_input_base,
        },
        setup=reaction_setup,
        conditions=reaction_conditions,
        workups=reaction_workups,
        outcomes=[outcome],
        provenance=ReactionProvenance(
            doi="10.1038/s41467-023-42446-5",
            record_created=RecordEvent(
                time=DateTime(value=datetime.datetime.now(tz=datetime.UTC).strftime("%Y-%m-%d %H:%M:%S %Z")),
                person=Person(username="qai", name="Qianxiang Ai", email="qai@mit.edu"))
        )
    )
    return reaction


if __name__ == '__main__':
    reactions = [create_reaction(r) for r in REACTION_RECORDS]
    dataset = Dataset(
        name="Dataset from `Ultra-high-throughput mapping of the chemical space of asymmetric catalysis enables accelerated reaction discovery`",
        description="Dataset from `Ultra-high-throughput mapping of the chemical space of asymmetric catalysis enables accelerated reaction discovery`",
        reactions=reactions,
    )
    write_message(dataset, "dataset.pbtxt")
