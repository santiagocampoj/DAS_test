// -----------------------------------------------------------------------------
// Demostración del efecto de ventana de la gauge length en DAS.
//
// Una onda acústica plana incide sobre la fibra con un ángulo theta. El DAS no
// mide la deformación en un punto, sino la diferencia finita del desplazamiento
// entre los dos extremos de un tramo de longitud Lg (la gauge length):
//
//     eps_das(z) = [ u(z + Lg/2) - u(z - Lg/2) ] / Lg
//
// Eso equivale a promediar la deformación verdadera sobre una ventana espacial
// rectangular de anchura Lg, lo que en el dominio del numero de onda axial se
// traduce en multiplicar por sinc(k_z * Lg / 2). Este programa lo comprueba.
//
// Compilar:  g++ -O2 -std=c++17 das_gauge.cpp -o das_gauge
// Ejecutar:  ./das_gauge                 (usa los parametros por defecto)
//            ./das_gauge <f> <theta> <Lg> <dz>   (sobrescribe sin recompilar)
//              f     = frecuencia [Hz]
//              theta = angulo de incidencia [grados]
//              Lg    = gauge length [m]
//              dz    = espaciado entre canales [m]
// -----------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <vector>

struct Params {
    double c = 1500.0;  // velocidad del sonido en agua [m/s]
    double f = 20.0;    // frecuencia de la senal [Hz]  (~ rorcual comun)
    double theta_deg = 45.0;    // angulo de incidencia [grados]
    double A = 1e-9;    // amplitud de desplazamiento [m]
    double Lg = 10.0;    // gauge length [m]
    double dz = 2.0;     // espaciado entre canales [m]
    double fiber_len = 2000.0;  // longitud de fibra [m]
    double noise = 0.0;     // ruido gaussiano (fraccion de la amplitud del strain); 0 = sin ruido
    long   seed = 42;      // semilla del ruido
};

static const double PI = 3.14159265358979323846;

// sinc no normalizada: sin(x)/x, con el limite sinc(0) = 1
static double sinc(double x) { return (x == 0.0) ? 1.0 : std::sin(x) / x; }

int main(int argc, char** argv) {
    Params p;
    if (argc > 1) p.f = std::atof(argv[1]);
    if (argc > 2) p.theta_deg = std::atof(argv[2]);
    if (argc > 3) p.Lg = std::atof(argv[3]);
    if (argc > 4) p.dz = std::atof(argv[4]);

    const double theta = p.theta_deg * PI / 180.0; // angulo de incidencia en radianes formula
    const double ka = 2.0 * PI * p.f / p.c;   // numero de onda acustico total: ka = 2 pi f / c
    const double kz = ka * std::sin(theta);   // numero de onda axial (proyectado sobre la fibra): kz = ka sin(theta)
    const double lambda_z = (kz == 0.0) ? INFINITY : 2.0 * PI / kz; // long. de onda axial aparente: lambda_z = 2 pi / kz

    // Campo, evaluado en un instante fijo t = 0 (una "traza" espacial del DAS).
    // Desplazamiento:      u(z)   = A sin(-kz z)
    // Strain verdadero:    eps(z) = du/dz = -A kz cos(-kz z)
    // Strain DAS:          diferencia finita de u sobre Lg, dividida por Lg: u(z+Lg/2) - u(z-Lg/2) / Lg
    auto u = [&](double z) { return p.A * std::sin(-kz * z); };
    auto eps_true = [&](double z) { return -p.A * kz * std::cos(-kz * z); };
    auto eps_das  = [&](double z) {
        return (u(z + p.Lg / 2.0) - u(z - p.Lg / 2.0)) / p.Lg;
    };

    std::mt19937 rng(p.seed); // generador de numeros aleatorios para el ruido
    const double strain_amp = p.A * kz; // amplitud del strain verdadero
    std::normal_distribution<double> gauss(0.0, p.noise * strain_amp);

    // ---- Salida 1: seccion espacial (datos "DAS" a lo largo de la fibra) ----
    {
        std::ofstream out("das_section.csv");
        out << "z_m,strain_verdadero,strain_das,strain_das_ruidoso\n";
        for (double z = 0.0; z <= p.fiber_len; z += p.dz) {
            double et = eps_true(z);
            double ed = eps_das(z);
            double en = ed + (p.noise > 0.0 ? gauss(rng) : 0.0);
            out << z << "," << et << "," << ed << "," << en << "\n";
        }
    }

    // ---- Salida 2: barrido de Lg/lambda_z frente a la sinc teorica ----
    // Para cada relacion r = Lg/lambda_z fijamos Lg = r*lambda_z y medimos, en un
    // punto donde el strain verdadero es maximo (z = 0), el cociente das/verdadero.
    // Debe coincidir con sinc(kz*Lg/2) incluyendo el signo de los lobulos laterales.
    {
        std::ofstream out("sweep.csv");
        out << "Lg_sobre_lambda,Lg_m,cociente_medido,sinc_teorico,error_abs\n";
        for (int i = 0; i <= 200; ++i) {

            double r  = i / 100.0; // relacion Lg/lambda_z de 0 a 2 
            double Lg = r * lambda_z;
            
            // ed es el strain medido por DAS en z=0, et es el strain verdadero en z=0
            double ed = (u(0.0 + Lg / 2.0) - u(0.0 - Lg / 2.0)) / Lg; // das en z=0
            double et = eps_true(0.0);                                 // verdadero en z=0
            
            double medido = (Lg == 0.0) ? 1.0 : ed / et;
            double teorico = sinc(kz * Lg / 2.0);
            
            out << r << "," << Lg << "," << medido << "," << teorico
                << "," << std::fabs(medido - teorico) << "\n";
        }
    }


    // ---- Resumen por pantalla ----
    double op_ratio = p.Lg / lambda_z; // punto de operacion actual
    double op_resp  = sinc(kz * p.Lg / 2.0); // respuesta del DAS en el punto de operacion actual
    const char* zona = (op_ratio < 0.2)  ? "BUENA (zona plana, sin distorsion)"
                     : (op_ratio < 0.44) ? "TRANSICION (atenuacion apreciable)"
                                         : "MALA (fuerte distorsion / cerca del cero)";


                                         
    std::printf("=== Parametros ===\n");
    std::printf("  c = %.1f m/s | f = %.2f Hz | theta = %.1f deg | Lg = %.2f m | dz = %.2f m\n",
                p.c, p.f, p.theta_deg, p.Lg, p.dz);
    std::printf("=== Derivados ===\n");
    std::printf("  longitud de onda en el medio  lambda   = %.2f m\n", p.c / p.f);
    std::printf("  longitud de onda axial        lambda_z = %.2f m\n", lambda_z);
    std::printf("  numero de onda axial          kz       = %.5f rad/m\n", kz);
    std::printf("=== Punto de operacion ===\n");
    std::printf("  Lg/lambda_z   = %.3f\n", op_ratio);
    std::printf("  respuesta sinc = %.4f  (pierdes ~%.1f%% de amplitud)\n",
                op_resp, 100.0 * (1.0 - std::fabs(op_resp)));
    std::printf("  zona: %s\n", zona);
    std::printf("=== Ficheros escritos ===\n");
    std::printf("  das_section.csv  (traza espacial: strain verdadero vs DAS)\n");
    std::printf("  sweep.csv        (barrido Lg/lambda_z: medido vs sinc)\n");
    return 0;
}