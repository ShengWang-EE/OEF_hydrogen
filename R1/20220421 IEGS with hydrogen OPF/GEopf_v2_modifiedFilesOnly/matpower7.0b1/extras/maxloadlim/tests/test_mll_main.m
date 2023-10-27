function success = test_mll_main(verbose, exit_on_fail)
%TEST_MLL_MAIN This function defines the list of tests to be run.
% The tests are run by using T_RUN_TESTS in MATPOWER.

%   MATPOWER
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Camille Hamon
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    verbose = 0;
end

tests = {};
tests{end+1} = 't_cpf_case9';
tests{end+1} = 't_cpf_case39';
tests{end+1} = 't_varGen_case9';
tests{end+1} = 't_varGen_case39';
t_run_tests( tests, verbose );
