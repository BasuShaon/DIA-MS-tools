
configfile: "driver.yaml"
import os; os.makedirs("output/logs", exist_ok=True) 

rule all:
    input:
        # Objects
        "output/shogoki.pkl", # carrier object
        "output/nigoki.pkl", # detectioncontroled carrier
        expand("output/sangoki_{bound}.pkl", bound = config["bound"]), # imputed carrier
        expand("output/yongoki_{bound}.pkl", bound = config["bound"]), # combat corrected carrier
        expand("output/gogoki_{bound}.pkl", bound = config["bound"]), # roll-downed carrier
        
        # Done files
        "output/logs/carrier.done",
        "output/logs/detectioncontrol.done",
        expand("output/logs/imputer_{bound}.done", bound = config["bound"]),
        expand("output/logs/combat_{bound}.done", bound = config["bound"]), 
        expand("output/logs/rolldown_{bound}.done", bound = config["bound"]) 

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

rule imputer: # Fan-out point
    input:
        script = "code/imputer.py",
        nigoki = "output/nigoki.pkl"
    params:
        bound = "{bound}"
    output:
        sangoki = "output/sangoki_{bound}.pkl",
        log_done = "output/logs/imputer_{bound}.done"
    log:
        "output/logs/imputer_{bound}.log"
    shell:
        """
        python3 {input.script} -b {params.bound} > {log} 2>&1
        touch {output.log_done}
        """

rule combat:
    input:
        script = "code/combat.py",
        sangoki = "output/sangoki_{bound}.pkl"
    output:
        yongoki = "output/yongoki_{bound}.pkl",
        log_done = "output/logs/combat_{bound}.done"
    log:
        "output/logs/combat_{bound}.log"
    shell:
        """
        python3 {input.script} -i {input.sangoki} > {log} 2>&1
        touch {output.log_done}
        """

rule rolldown:
    input:
        script = "code/rolldown.py",
        yongoki = "output/yongoki_{bound}.pkl"
    output:
        gogoki = "output/gogoki_{bound}.pkl",
        log_done = "output/logs/rolldown_{bound}.done"
    log:
        "output/logs/rolldown_{bound}.log"
    shell:
        """
        python3 {input.script} -i {input.yongoki} > {log} 2>&1
        touch {output.log_done}
        """
