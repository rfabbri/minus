-- Utilities to help iterface with Minus
-- To be included ih other files

-- Version of cCode that outputs to file other than stdout ---------------------------------------
cCode (String, List,List) := (outfile, outputs,inputs) -> (
    f := openOut outfile;
    h := newPrintTable " = ";
    scan(inputs, g->printName(g,h));
    scan(outputs, g->printName(g,h));
    f << "{" << endl;
    scan(h#"#vars", i-> f << ("  C<F> &X"|i|" = x["|i|"];") << endl);
    scan(h#"#lines", i-> f << ("  const C<F> "| h#i | ";") << endl);
    scan(#outputs, i-> f << ("  y["|i|"] = "|printName(outputs#i,h)|";") << endl); 
    f << "}";
    close f;
    )

cCode (String, GateMatrix,GateMatrix) := (outfile, M,I) -> cCode(outfile,
flatten entries M, flatten entries I)
