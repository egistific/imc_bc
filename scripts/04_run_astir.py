import argparse
import numpy as np
import pandas as pd

# astir 0.1.5 uses np.NaN / np.Inf which have been removed in numpy 2.4
np.NaN = np.nan
np.Inf = np.inf

from astir.data import from_csv_yaml


p = argparse.ArgumentParser(description="Run Astir cell phenotyping")
p.add_argument("--exprs", required=True, type=str, help="expression csv")
p.add_argument("--design", required=False, default=None, help="design csv")
p.add_argument("--markers", required=True, type=str, help="markers yaml")
p.add_argument("--output", required=True, type=str, help="output csv")
p.add_argument("--seed", default=42, type=int, help="seed")
p.add_argument("--max_epochs", default=50, type=int, help="maximum number of epochs")
p.add_argument("--batch_size", default=1024, type=int, help="batch size")
p.add_argument("--learning_rate", default=1e-3, type=float, help="learning_rate")
p.add_argument("--n_init", default=5, type=int, help="number of models initially trained")
p.add_argument("--n_init_epochs", default=5, type=int, help="number of epochs initial models trained for")
p.add_argument("--delta_loss", default=1e-3, type=float, help="delta loss")
p.add_argument("--threshold", default=0.7, type=float, help="threshold")

args = p.parse_args()

if args.design is not None:
    # Check that cells are in the same order in both files
    ids_expr   = pd.read_csv(args.exprs,  usecols=[0]).iloc[:, 0]
    ids_design = pd.read_csv(args.design, usecols=[0]).iloc[:, 0]
    assert ids_expr.equals(ids_design), "Expression and design rows do not match!"
    ast = from_csv_yaml(args.exprs, marker_yaml=args.markers,
                        design_csv=args.design, random_seed=args.seed)
else:
    ast = from_csv_yaml(args.exprs, marker_yaml=args.markers, random_seed=args.seed)
    print("No batch effects correction!")
print(ast)

# Train model
ast.fit_type(max_epochs=args.max_epochs, batch_size=args.batch_size, learning_rate=args.learning_rate,
            n_init=args.n_init, n_init_epochs=args.n_init_epochs, delta_loss=args.delta_loss)

losses = ast.get_type_losses()
print(losses[-5:])

print(ast.diagnostics_celltype())

prob = ast.get_celltype_probabilities()
labels = ast.get_celltypes(threshold=args.threshold).iloc[:, 0]
print(labels.value_counts(normalize=True).round(3))

prob = ast.get_celltype_probabilities()
print(prob.max(axis=1).describe()) # take each cell's highest probability and get the summary

out = prob.copy()
out["astir_label"] = labels
out.to_csv(args.output)
