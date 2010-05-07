
import enrichment
import sys

def ratio (b, n, B, N):
	return (float(b) / n) / (float(B) / N)

def print_contingency_table (b, n, B, N, stream = sys.stdout):
	print >>stream, "  %5s (b)    %5s      | %5s (n)" % (b, n - b, n)
	print >>stream, "  %5s        %5s      | %5s" % (B - b, N - B - n + b, N - n)
	print >>stream, "  ------------------------+-------"
	print >>stream, "  %5s (B)    %5s      | %5s (N)" % (B, N - B, N)

	print >>stream, "\n              Ratio: %.6f" % ratio(b, n, B, N)

	left, right, two_tailed = enrichment.fisher_exact_test(b, n, B, N)

	print >>stream, "       Left p-value: %.6g" % left
	print >>stream, "      Right p-value: %.6g" % right
	print >>stream, " Two-tailed p-value: %.6g" % two_tailed
