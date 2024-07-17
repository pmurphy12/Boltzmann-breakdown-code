function KLD = KLtest_continuous(P,Q,dx)

    KLD = sum( P .* log2( P./Q ) * dx , 'omitnan');

end