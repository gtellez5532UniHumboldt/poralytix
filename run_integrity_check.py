#!/usr/bin/env python3
"""
Script de Verificación de Integridad Integral del Proyecto Refactorizado.
Ejecuta diagnósticos de entorno, estructura de paquetes, pruebas unitarias
y validación de datasets experimentales de prueba.

Uso:
    py -3.11 run_integrity_check.py
"""
import sys
import os
import time
import unittest
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))


def print_banner(text):
    print("\n" + "=" * 68)
    print(f"  {text}")
    print("=" * 68)


def check_environment():
    print_banner("1. DIAGNÓSTICO DE ENTORNO PYTHON Y DEPENDENCIAS")
    print(f"  • Versión de Python:    {sys.version.split()[0]}")
    print(f"  • Ejecutable:           {sys.executable}")
    print(f"  • Plataforma:           {sys.platform}")

    packages = ['numpy', 'scipy', 'pandas', 'matplotlib', 'reportlab', 'tkinter']
    missing = []
    for pkg in packages:
        try:
            mod = __import__(pkg)
            ver = getattr(mod, '__version__', 'N/A')
            print(f"  [OK] {pkg:<14} v{ver}")
        except Exception as e:
            print(f"  [ERROR] {pkg:<14} -> NO DISPONIBLE ({e})")
            missing.append(pkg)

    if missing:
        print(f"\n[ALERTA] Faltan dependencias: {missing}")
        return False
    return True


def check_test_datasets():
    print_banner("2. VERIFICACIÓN DE DATASETS EXPERIMENTALES DE PRUEBA")
    root_dir = BASE_DIR.parent
    datasets = ['813', '820', 'sat1', 'CZS5c-']
    all_ok = True

    for ds in datasets:
        p = root_dir / ds
        if p.exists() and p.is_dir():
            files = list(p.glob("**/*.fcs")) + list(p.glob("**/*.FCS"))
            total_bytes = sum(f.stat().st_size for f in files)
            print(f"  [OK] Dataset {ds:<10} -> {len(files):>2} archivos .fcs ({total_bytes / (1024*1024):.1f} MB)")
        else:
            print(f"  [ERROR] Dataset {ds} no encontrado en {root_dir}")
            all_ok = False
    return all_ok


def run_unit_tests():
    print_banner("3. EJECUCIÓN DE SUITE DE PRUEBAS AUTOMATIZADAS")
    tests_dir = BASE_DIR / "tests"
    loader = unittest.TestLoader()
    suite = loader.discover(start_dir=str(tests_dir), pattern="test_*.py")

    start_time = time.time()
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    elapsed = time.time() - start_time

    print(f"\nResumen de Pruebas:")
    print(f"  • Total de tests ejecutados: {result.testsRun}")
    print(f"  • Errores:                   {len(result.errors)}")
    print(f"  • Fallos:                    {len(result.failures)}")
    print(f"  • Tiempo transcurrido:       {elapsed:.3f} s")

    return result.wasSuccessful()


def main():
    print_banner("ESTACIÓN DE TRABAJO DE CITOMETRÍA — CONTROL DE INTEGRIDAD")
    env_ok = check_environment()
    data_ok = check_test_datasets()
    tests_ok = run_unit_tests()

    print_banner("4. DICTAMEN FINAL DE INTEGRIDAD")
    if env_ok and data_ok and tests_ok:
        print("  >>> ESTADO GENERAL: [APROBADO] - SISTEMA 100% OPERATIVO <<<")
        print("  • Todos los módulos I/O, matemáticos, GUI y exportación validados.")
        print("  • Datasets de prueba 813, 820, sat1 y CZS5c- íntegros.")
        print("  • Heurística de metadatos para CZS5c- y subcarpetas confirmada.")
        sys.exit(0)
    else:
        print("  >>> ESTADO GENERAL: [REVISIÓN REQUERIDA] <<<")
        sys.exit(1)


if __name__ == "__main__":
    main()
