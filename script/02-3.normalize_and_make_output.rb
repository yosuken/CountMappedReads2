
fcnt, fcov, frlt = ARGV

# OBV_N00003_1 OBV_N00003 3 278 + 276 0
header       = %w|feature sequence start end strand length count coverage FPKM(RPKM) TPM|
t_cnt        = 0
t_cnt_per_bp = 0
gid2out      = {}

IO.readlines(fcnt)[2..-1].each{ |l|
	gid, sid, start, stop, strand, len, cnt = l.chomp.split("\t")
	raise if gid2out[gid]

	cnt    = cnt.to_i
	len    = len.to_i
	t_cnt += cnt                   ## for FPKM(RPKM)
	cnt_per_bp    = cnt.to_f / len ## for TPM
	t_cnt_per_bp += cnt_per_bp     ## for TPM

	gid2out[gid]  = [gid, sid, start, stop, strand, len, cnt, 0, cnt, cnt_per_bp]
}

IO.readlines(fcov).each{ |l|
	sid, start, stop, gid, score, strand, sumcov = l.chomp.split("\t")
	gid2out[gid][7] += sumcov.to_i
}

gid2out.each{ |gid, a|
	len    = a[5]
	sumcov = a[7]
	cov    = sumcov.to_f / len
	fpkm   = a[8] * 1_000.0 / len * 1_000_000 / t_cnt
	tpm    = a[9] * 1_000_000.0 / t_cnt_per_bp
	gid2out[gid] = a[0..6] + [cov, fpkm, tpm].map{ |i| "%.3f" % i }
}

open(frlt, "w"){ |fout|
	fout.puts header*"\t"
	gid2out.each{ |gid, a|
		fout.puts a*"\t"
	}
}
