fn polynomial() {
    let x_nodes: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y_nodes: Vec<f64> = x_nodes
        .iter()
        .map(|&x| 3.0 * x.powi(8) - 2.0 * x.powi(5) + 3.0 * x.powi(2) + 1.0)
        .collect();

    let mut matrix = vec![vec![0.0; 5]; 5];
    let mut rhs = vec![0.0; 5];

    for (i, &x) in x_nodes.iter().enumerate() {
        for j in 0..5 {
            matrix[i][j] = x.powi(j as i32);
        }
        rhs[i] = y_nodes[i];
    }

    let coefficients = solve_system(matrix, rhs);
    println!("Polynomial coefficients: {:?}", coefficients);
}

// Розв'язання системи методом Гаусса
fn solve_system(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = matrix.len();
    for i in 0..n {
        // Прямий хід
        let pivot = matrix[i][i];
        for j in i..n {
            matrix[i][j] /= pivot;
        }
        rhs[i] /= pivot;

        for k in (i + 1)..n {
            let factor = matrix[k][i];
            for j in i..n {
                matrix[k][j] -= factor * matrix[i][j];
            }
            rhs[k] -= factor * rhs[i];
        }
    }

    // Зворотній хід
    let mut solution = vec![0.0; n];
    for i in (0..n).rev() {
        solution[i] = rhs[i];
        for j in (i + 1)..n {
            solution[i] -= matrix[i][j] * solution[j];
        }
    }

    solution
}

fn cubic_spline() {
    let x_nodes: Vec<f64> = vec![0.0, 2.0, 4.0];
    let y_nodes: Vec<f64> = x_nodes
        .iter()
        .map(|&x| 3.0 * x.powi(8) - 2.0 * x.powi(5) + 3.0 * x.powi(2) + 1.0)
        .collect();

    let derivatives = vec![
        derivative(0.0),
        derivative(4.0),
    ];

    let spline_coefficients = calculate_cubic_spline(&x_nodes, &y_nodes, &derivatives);
    println!("Cubic spline coefficients: {:#?}", spline_coefficients);
}

fn calculate_cubic_spline(x: &[f64], y: &[f64], derivatives: &[f64]) -> Vec<Vec<f64>> {
    let n = x.len() - 1;
    let mut coeffs = vec![vec![0.0; 4]; n];

    for i in 0..n {
        let h = x[i + 1] - x[i];
        let a = y[i];
        let b = derivatives[i];
        let c = if i < n - 1 {
            (3.0 * (y[i + 1] - y[i]) / h - 2.0 * derivatives[i] - derivatives[i + 1]) / h
        } else {
            0.0 // Assuming zero or a different boundary condition for the last interval
        };

        let d = if i < n - 1 {
            (2.0 * (y[i] - y[i + 1]) / h + derivatives[i] + derivatives[i + 1]) / h.powi(2)
        } else {
            0.0 // Assuming zero or a different boundary condition for the last interval
        };

        coeffs[i] = vec![a, b, c, d];
    }

    coeffs
}

fn derivative(x: f64) -> f64 {
    24.0 * x.powi(7) - 10.0 * x.powi(4) + 6.0 * x
}

fn linear_spline() {
    let start = 0.0;
    let end = 4.0;
    let step = 0.5;

    let x_nodes: Vec<f64> = (0..)
        .map(|i| start + i as f64 * step)
        .take_while(|&x| x <= end)
        .collect();

    let y_nodes: Vec<f64> = x_nodes
        .iter()
        .map(|&x| 3.0 * x.powi(8) - 2.0 * x.powi(5) + 3.0 * x.powi(2) + 1.0)
        .collect();

    for i in 0..(x_nodes.len() - 1) {
        let x1 = x_nodes[i];
        let x2 = x_nodes[i + 1];
        let y1 = y_nodes[i];
        let y2 = y_nodes[i + 1];
        let slope = (y2 - y1) / (x2 - x1);

        println!(
            "Linear spline [{}, {}]: S(x) = {:.4} + {:.4}(x - {:.4})",
            x1, x2, y1, slope, x1
        );
    }
}


fn main() {
    polynomial();
    cubic_spline();
    linear_spline();
}
