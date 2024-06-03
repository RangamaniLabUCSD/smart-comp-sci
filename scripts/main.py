from typing import NamedTuple, Callable
import arguments
import runner

class Command(NamedTuple):
    msg: str
    script: Callable


commands = {
    "convert-notebooks": Command(
        msg="Convert notebook",
        script=runner.convert_notebooks
    ),
    "dendritic-spine-preprocess": Command(
        msg="Run preprocess dendritic spine mesh",
        script=runner.preprocess_spine_mesh
    ),
    "dendritic-spine": Command(
        msg="Run dendritic spine example",
        script=runner.dendritic_spine_example
    ),
    "dendritic-spine-postprocess": Command(
        msg="Postprocess dendritic spine example",
        script=runner.dendritic_spine_postprocess
    ),
    "cru-preprocess": Command(
        msg="Run preprocess cru mesh",
        script=runner.preprocess_cru_mesh
    ),
    "cru": Command(
        msg="Run cru example",
        script=runner.cru_example
    ),
    "mechanotransduction-preprocess": Command(
        msg="Run preprocess mechanotransduction mesh",
        script=runner.pre_preprocess_mech_mesh
    ),
    "mechanotransduction": Command(
        msg="Run mechanotransduction example",
        script=runner.mechanotransduction_example
    ),
    "mechanotransduction-postprocess": Command(
        msg="Postprocess mechanotransduction example",
        script=runner.postprocess_mechanotransduction
    ),
    "mito-preprocess": Command(
        msg="Run preprocess mito mesh",
        script=runner.preprocess_mito_mesh
    ),
    "mito": Command(
        msg="Run mito example",
        script=runner.mito_example
    ),
    "phosphorylation-preprocess": Command(
        msg="Run preprocess phosphorylation mesh",
        script=runner.preprocess_phosphorylation_mesh
    ),
    "phosphorylation": Command(
        msg="Run phosphorylation example",
        script=runner.phosphorylation_example
    ),
    "phosphorylation-postprocess": Command(
        msg="Postprocess phosphorylation example",
        script=runner.postprocess_phosphorylation
    ),
}

def main():
    parser = arguments.setup_parser()
    args = vars(parser.parse_args())

    try:
        cmd = commands[args["command"]]
        print(cmd.msg)
        cmd.script(**args)
    except KeyError:
        raise ValueError(f"Invalid command: {args['command']}")


if __name__ == "__main__":
    raise SystemExit(main())
