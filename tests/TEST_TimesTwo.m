% setup

%% Test 1: input = 2
in = 2;
assert(TimesTwo(in) == 4)

%% Test 2: input = -1
in = -1;
assert(TimesTwo(in) == -2)

%% Test 3: now passes!
in = 0;
assert(TimesTwo(in) == 2 * in) 