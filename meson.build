project('checkDISgrad', 'cpp')
nnpdf = dependency('nnpdf')
lhapdf = dependency('lhapdf')
e = executable('checkDISgrad', 'main.cc', 'toypdf/toy.cc', dependencies: [nnpdf, lhapdf])
test('checkDISgrad', e)
