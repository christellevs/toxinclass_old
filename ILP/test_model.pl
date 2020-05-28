:- use_module('source/gilps').

:- show_settings.
:- set(verbose, 3).
proteins:-
  read_problem('datasets/proteins/proteins'),
  evaluate_theory('models/theory_proteins.pl', 'datasets/proteins/examples_test').

:-proteins.