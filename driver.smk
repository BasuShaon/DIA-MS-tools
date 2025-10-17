configfile: "params.yaml"
import os; os.makedirs("output/logs", exist_ok=True) 

rule all:
    input:
        "output/shogoki.pkl",
        "output/nigoki.pkl",
        "output/sangoki.pkl",
        "output/yongoki.pkl",
        "output/gogoki.pkl",
        "output/logs/carrier.done",
        "output/logs/detectioncontrol.done",
        "output/logs/imputer.done",
        "output/logs/combat.done",
        "output/logs/rollup.done"

rule init_carrier:
    input:
        script = "code/carrier.py"
    output:
        shogoki = "output/shogoki.pkl",
        log_done = "output/logs/carrier.done"
    log:
        "output/logs/carrier.log"
    params:
        data_df = config["data_df"],
        meta_df = config["meta_df"],
        proj = config["proj"]
    shell:
        """
        python3 {input.script} -d {params.data_df} -m {params.meta_df} -p {params.proj} > {log} 2>&1
        touch {output.log_done}
        """

rule detectioncontrol:
    input:
        script = "code/detectioncontrol.py",
        shogoki = "output/shogoki.pkl"
    output:
        nigoki = "output/nigoki.pkl",
        log_done = "output/logs/detectioncontrol.done"
    log:
        "output/logs/detectioncontrol.log"
    params:
        alpha = config["alpha"]
    shell:
        """
        python3 {input.script} -a {params.alpha} > {log} 2>&1
        touch {output.log_done}
        """

rule imputer:
    input:
        script = "code/imputer.py",
        nigoki = "output/nigoki.pkl"
    output:
        sangoki = "output/sangoki.pkl",
        log_done = "output/logs/imputer.done"
    log:
        "output/logs/imputer.log"
    params:
        bound = config["bound"]
    shell:
        """
        python3 {input.script} -b {params.bound} > {log} 2>&1
        touch {output.log_done}
        """

rule combat:
    input:
        script = "code/combat.py",
        sangoki = "output/sangoki.pkl"
    output:
        yongoki = "output/yongoki.pkl",
        log_done = "output/logs/combat.done"
    log:
        "output/logs/combat.log"
    shell:
        """
        python3 {input.script} > {log} 2>&1
        touch {output.log_done}
        """

rule rollup:
    input:
        script = "code/rollup.py",
        yongoki = "output/yongoki.pkl"
    output:
        gogoki = "output/gogoki.pkl",
        log_done = "output/logs/rollup.done"
    log:
        "output/logs/rollup.log"
    shell:
        """
        python3 {input.script} > {log} 2>&1
        touch {output.log_done}
        """
