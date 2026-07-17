import gurobipy as gp
from gurobipy import GRB
import numpy as np

def benders_decomposition(A, B, b, c, f, x_bounds, y_bounds, max_iters=50, tol=1e-6):

    m, nx = A.shape
    _, ny = B.shape

    # Initialize Master Problem
    master = gp.Model("Benders_Master")
    x = master.addVars(nx, lb=[lb for lb, ub in x_bounds],
                          ub=[ub for lb, ub in x_bounds],
                          name="x", vtype=GRB.CONTINUOUS)
    theta = master.addVar(name="theta", vtype=GRB.CONTINUOUS)
    master.setObjective(gp.quicksum(c[i]*x[i] for i in range(nx)) + theta, GRB.MINIMIZE)

    cuts = []

    # Benders Loop
    for it in range(max_iters):
        master.optimize()

        if master.Status != GRB.OPTIMAL:
            print("Master problem not optimal!")
            break

        x_val = np.array([x[i].X for i in range(nx)])

        # Subproblem: fixed x
        sub = gp.Model("Benders_Subproblem")
        y = sub.addVars(ny, lb=[lb for lb, ub in y_bounds],
                           ub=[ub for lb, ub in y_bounds],
                           name="y", vtype=GRB.CONTINUOUS)

        # Constraints: A x + B y >= b => B y >= b - A x
        rhs = b - A @ x_val
        constraints = [sub.addConstr(gp.quicksum(B[i, j] * y[j] for j in range(ny)) >= rhs[i])
                       for i in range(m)]

        sub.setObjective(gp.quicksum(f[j] * y[j] for j in range(ny)), GRB.MINIMIZE)
        sub.setParam("OutputFlag", 0)
        sub.optimize()

        if sub.Status == GRB.OPTIMAL:
            y_val = np.array([y[j].X for j in range(ny)])
            sub_obj = f @ y_val
            theta_val = theta.X

            if theta_val + tol < sub_obj:
                # Add optimality cut
                duals = np.array([c.Pi for c in constraints])
                cut_expr = gp.LinExpr()
                cut_expr += sum(duals[i] * (b[i] - gp.quicksum(A[i, k] * x[k] for k in range(nx)))
                                for i in range(m))
                master.addConstr(theta >= cut_expr)
                cuts.append("optimality")
            else:
                print(f"Converged in {it+1} iterations.")
                return x_val, y_val, c @ x_val + f @ y_val
        elif sub.Status == GRB.INFEASIBLE:
            # Feasibility cut via dual Farkas
            sub.computeIIS()
            sub.write("infeasible.ilp")
            sub.Params.InfUnbdInfo = 1
            sub.optimize()
            duals = np.array([c.FarkasDual for c in constraints])
            cut_expr = sum(duals[i] * (b[i] - gp.quicksum(A[i, k] * x[k] for k in range(nx)))
                           for i in range(m))
            master.addConstr(cut_expr <= 0)
            cuts.append("feasibility")
        else:
            print("Subproblem not solved properly.")
            break

    print("Benders did not converge.")
    return None, None, None


# Example problem data
A = np.array([[1, 2], [3, 4]])
B = np.array([[1, 0], [0, 1]])
b = np.array([5, 11])
c = np.array([1.0, 2.0])
f = np.array([3.0, 1.0])
x_bounds = [(0, GRB.INFINITY), (0, GRB.INFINITY)]
y_bounds = [(0, GRB.INFINITY), (0, GRB.INFINITY)]

x_sol, y_sol, obj = benders_decomposition(A, B, b, c, f, x_bounds, y_bounds)
print(f"Optimal x: {x_sol}, y: {y_sol}, Objective: {obj}")
