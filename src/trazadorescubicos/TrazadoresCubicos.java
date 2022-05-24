package trazadorescubicos;

public class TrazadoresCubicos {

    public static void main(String[] args) {
        double[] x = {0, 1, 2, 3, 4};
        double[] y = {4, 2, -2, -2, 8};
        
        int n = x.length;
        
        //se crea la matriz de los trazadores
        double[][] a = new double[n - 2][n - 2];
        double[] b = new double[n - 2];
        
        for (int i = 1; i < n - 1; i++) {
            if (i > 1) {
                a[i - 1][i - 2] = x[i] - x[i - 1];
            }
            a[i - 1][i - 1] = 2 * (x[i + 1] - x[i - 1]);
            if (i < n - 2) {
                a[i - 1][i] = x[i + 1] - x[i];
            }
            b[i - 1] = (6 / (x[i + 1] - x[i])) * (y[i + 1] - y[i]) + (6 / (x[i] - x[i - 1]) * (y[i - 1] - y[i]));
        }
                
        double[] rtaAux = gaussJordan(a, b);
        double[] f2 = new double[n];
        System.arraycopy(rtaAux, 0, f2, 1, rtaAux.length);
        
        for (int i = 1; i < n; i++) {
            double t1 = f2[i - 1] / (6 * (x[i] - x[i - 1]));
            double t2 = f2[i] / (6 * (x[i] - x[i - 1]));
            double t3 = y[i - 1] / (x[i] - x[i - 1]) - f2[i - 1] * (x[i] - x[i - 1]) / 6;
            double t4 = y[i] / (x[i] - x[i - 1]) - f2[i] * (x[i] - x[i - 1]) / 6;
            
            double[] arrCoef = new double[4];
            
            //Se calculan los coeficientes del polinomio
            arrCoef[0] = t1 * Math.pow(x[i], 3) - t2 * Math.pow(x[i - 1], 3) + t3 * x[i] - t4 * x[i - 1];
            arrCoef[1] = -t1 * 3 * Math.pow(x[i], 2) + t2 * 3 * Math.pow(x[i - 1], 2) - t3 + t4;
            arrCoef[2] = t1 * 3 * x[i] - t2 * 3 * x[i - 1];
            arrCoef[3] = -t1 + t2;
            
            System.out.print("f(x) = ");
            for (int j = 0; j < 4; j++) {
                if (arrCoef[j] != 0) {
                    if (j > 0) {
                        System.out.print("+ ");
                    }
                    System.out.print(arrCoef[j]);
                    switch (j) {
                        case 0:
                            System.out.print(" ");
                            break;
                        case 1:
                            System.out.print("x ");
                            break;
                        default:
                            System.out.print("x^" + j + " ");
                            break;
                    }
                }
            }
            System.out.println("{x>=" + x[i - 1] + "}{x<" + x[i] + "}");
        }
    }
    
    private static double[] gaussJordan(double[][] a, double[] b) {
        double[][] aAux = duplicarArreglo(a);
        double[] bAux = duplicarArreglo(b);
        
        int n = bAux.length;
        
        //Se construye la matriz triangular superior
        for (int i = 0; i < n; i++) {
            //Pivoteo
            if (aAux[i][i] == 0) {
                for (int k = i + 1; k < n; k++) {
                    if (aAux[k][i] != 0) {
                        double[] filaAux = aAux[i];
                        aAux[i] = aAux[k];
                        aAux[k] = filaAux;
                        
                        double valoAux = bAux[i];
                        bAux[i] = bAux[k];
                        bAux[k] = valoAux;
                        break;
                    }
                }
            }
            
            //Escalonamiento
            double valorAux = aAux[i][i];
            for (int j = i; j < n; j++) {
                aAux[i][j] /= valorAux;
            }
            bAux[i] /= valorAux;
            
            //Reducción
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    double fact = aAux[j][i] / aAux[i][i];

                    for (int k = 0; k < n; k++) {
                        aAux[j][k] -= (aAux[i][k] * fact);
                    }
                    
                    bAux[j] -= (bAux[i] * fact);
                }
            }
        }
        
        return bAux;
    }
    
    private static double[][] duplicarArreglo(double[][] m) {
        double[][] duplicado = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            System.arraycopy(m[i], 0, duplicado[i], 0, m[i].length);
        }
        
        return duplicado;
    }
    
    private static double[] duplicarArreglo(double[] v) {
        double[] duplicado = new double[v.length];
        System.arraycopy(v, 0, duplicado, 0, v.length);
        
        return duplicado;
    }
    
}