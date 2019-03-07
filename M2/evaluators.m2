-- code for generating various evaluators 

-- HxHt
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})

-- HxH
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})
