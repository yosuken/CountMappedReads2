###
### CountMappedReads2.rake: a tool for NGS read mapping, count reads and calculate normalized values
###
### Copyright: 2018 (C) Yosuke Nishimura (ynishimura@aori.u-tokyo.ac.jp)
### Licence: MIT license
###


# {{{ procedures
WriteBatch  = lambda do |outs, jdir, t|
	outs.each_slice(10000).with_index(1){ |ls, idx|
		open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.#{idx}", "w"){ |fjob|
			fjob.puts ls
		}
	}
end

RunBatch    = lambda do |jdir, queue, nthreads, mem, wtime, ncpus|
	# [TODO] queue validation
	Dir["#{jdir}/*"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin|
		if queue != ""
			raise("`--queue #{queue}': invalid queue") unless %w|JP1 JP4 JP10 cdb|.include?(queue)
			sh "qsubarraywww -q #{queue} -l ncpus=#{nthreads} -l mem=#{mem}gb -l walltime=#{wtime} #{fin}"
		elsif ncpus != ""
			raise("`--ncpus #{ncpus}': not an integer") if ncpus !~ /^\d+$/
			sh "parallel --jobs #{ncpus} <#{fin}"
		else
			sh "sh #{fin}"
		end
	}
end

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "\e[1;32m===== #{Time.now}\e[0m"
	puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
	puts ""
	$stdout.flush
end

CheckVersion = lambda do |commands|
	commands.each{ |command|
		str = case command
					when "ruby"
						%|LANG=C ruby --version 2>&1|
					when "makeblastdb"
						%|LANG=C makeblastdb -version 2>&1|
					when "tblastx"
						%|LANG=C tblastx -version 2>&1|
					when "R"
						%{LANG=C R --version 2>&1 |head -n 3}
					when "ape"
						%|LANG=C R --quiet --no-save --no-restore -e "packageVersion('ape')" 2>&1|
					when "phangorn"
						%|LANG=C R --quiet --no-save --no-restore -e "packageVersion('phangorn')" 2>&1|
					when "DESeq2"
						%|LANG=C R --quiet --no-save --no-restore -e "packageVersion('DESeq2')" 2>&1|
					when "parallel"
						%{LANG=C parallel --version 2>&1 |head -n 1}
					when "bowtie2"
						%{LANG=C bowtie2 --help 2>&1 |head -n 1}
					when "samtools"
						%{LANG=C samtools --help 2>&1 |sed -ne '3p'}
					when "featureCounts"
						%{LANG=C featureCounts 2>&1 |sed -ne '2p'}
					end
		puts ""
		puts "\e[1;32m===== check version: #{command}\e[0m"
		puts ""
		puts "$ #{str}"
		### run
		puts `#{str}`
		### flush
		$stdout.flush
	}
end
# }}} procedures


# {{{ default (run all tasks)
task :default do
	### define shared tasks
	tasks = %w|
	01-1.bowtie2_build 01-2.bowtie2 01-3.samtools 01-4.gene_file_convert 01-5.make_sequence.gtf
	02-1.featureCounts 02-2.bedcov 02-3.normalize_and_make_output
	|

	### constants
	Odir     = ENV["dir"]     ## output directory
	Fin      = ENV["fin"]     ## input fasta file

	## featureCounts
	Gene     = ENV["gene"]    ## input gene annotation file
	Feature  = ENV["feature"] ## raise "--feature sequence is not allowed" if Feature == "sequence"
	Attr     = ENV["attr"] 
	if Gene == "" ## when gene annotaion is not given
		Types    = ["sequence"]
		Attrs    = ["sequence_id"]
		Format   = ""
	else
		Types    = [Feature, "sequence"]
		Attrs    = [Attr,    "sequence_id"]
		if ENV["format"].size > 0
			Format = ENV["format"]
		else
			Format = Gene.split(".")[-1].upcase ## given or inferred from extention of Gene
			## Format validation
			unless %w|GTF GTF2 GFF3 BED BED6|.include?(Format)
				raise "--format can not be inferred from .#{Gene.split(".")[-1]}"
			end
		end
	end

	## Bowtie2
	Pe1      = ENV["pe1"] ## paired-end reads 1
	Pe2      = ENV["pe2"] ## paired-end reads 2
	Up       = ENV["up"]  ## single-end reads
	ScoreMin = ENV["score_min"]
	# end-to-end, ~80% idt: "L,0,-1.2"
	# end-to-end, ~85% idt: "L,0,-0.9"
	# end-to-end, ~90% idt: "L,0,-0.6"
	# end-to-end, ~95% idt: "L,0,-0.3"

	## dir/file
	Idir     = "#{Odir}/sequence"     
	Fa       = "#{Idir}/#{File.basename(Fin)}"
	Mdir     = "#{Odir}/mapping"
	Adir     = "#{Odir}/annotation"
	Cdir     = "#{Odir}/count"
	Sam      = "#{Mdir}/out.sam"

	Threads  = ENV["threads"]||"1"
	Mem      = ENV["mem"]||"1G"

	### check version
	commands  = %w|bowtie2 samtools featureCounts ruby|
	# commands += %w|parallel| if Ncpus != ""
	CheckVersion.call(commands)

	### run
	NumStep  = tasks.size
	tasks.each.with_index(1){ |task, idx|
		Rake::Task[task].invoke(idx)
	}
end
# }}} default (run all tasks)


# {{{ tasks 01
desc "01-1.bowtie2_build"
task "01-1.bowtie2_build", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	mkdir_p Idir

	sh "cp #{Fin} #{Fa}"
	sh "bowtie2-build --threads #{Threads} #{Fa} #{Fa}"
end
desc "01-2.bowtie2"
task "01-2.bowtie2", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	mkdir_p Mdir

	input     = ""
	input    += "-1 #{Pe1} " if Pe1 != ""
	input    += "-2 #{Pe2} " if Pe2 != ""
	input    += "-U #{Up} "  if Up  != ""

	sh "bowtie2 --end-to-end --score-min #{ScoreMin} --threads #{Threads} -x #{Fa} #{input} -S #{Mdir}/out.sam 2>#{Mdir}/bowtie2.log"
end
desc "01-3.samtools"
task "01-3.samtools", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

	bam       = "#{Mdir}/out.bam.tmp"
	sorted    = "#{Mdir}/out.bam"

	str =  ["samtools view --threads #{Threads} -bS #{Sam} -o #{bam} 2>#{bam}.makelog",
				  "rm #{Sam}",
					"samtools sort --threads #{Threads} -m #{Mem} -O bam -T #{Mdir} #{bam} -o #{sorted} 2>#{sorted}.makelog",
				  "rm #{bam}",
					"samtools index #{sorted} 2>#{sorted}.indexlog",
	        ].join(" && ")

	sh "#{str}"
end
desc "01-4.gene_file_convert"
task "01-4.gene_file_convert", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	mkdir_p Adir

	if Gene == "" ## when gene annotaion is not given
	else
		case Format
		when "GTF", "GTF2"
			# <sid> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
			# [example of attributes] gene_id "001"; transcript_id "001.1";
			open("#{Adir}/#{Feature}.gtf", "w"){ |fout| fout.puts IO.read(Gene) }
			outs = []
			IO.readlines(Gene).each{ |l|
				next if l =~ /^#/
				next if l =~ /^$/
				sid, source, feature, start, stop, score, strand, frame, attrs = l.chomp.split("\t")
				name = attrs[/#{Attr} "([^"]+)"/, 1]
				next unless name
				outs << [sid, start.to_i-1, stop, name, score, strand]
			}
			open("#{Adir}/#{Feature}.bed", "w"){ |fout| 
				outs.each{ |a| fout.puts a*"\t" }
			}
		when "GFF3"
			# <sid> <source> <type(feature)> <start> <end> <score> <strand> <phase(frame)> [attributes]
			open("#{Adir}/#{Feature}.gff3", "w"){ |fout| fout.puts IO.read(Gene) }
			outs = []
			IO.readlines(Gene).each{ |l|
				next if l =~ /^#/
				next if l =~ /^$/
				sid, source, feature, start, stop, score, strand, frame, attrs = l.chomp.split("\t")
				name = attrs[/#{Attr}=([^;]+)/, 1]
				next unless name
				outs << [sid, start.to_i-1, stop, name, score, strand]
			}
			open("#{Adir}/#{Feature}.bed", "w"){ |fout| 
				outs.each{ |a| fout.puts a*"\t" }
			}
		when "GFF2"
			raise("Currently, GFF2 is not supported.")
			# GFF2 (MetaGeneMark style (gene_id=XXX))
			# <id>  GeneMark.hmm    CDS     90      560     18.620995       +       0       gene_id=1
		when "BED", "BED6"
			# <sid> <0-based-start> <end> <name> <score> <strand>
			outs1 = []
			outs2 = []
			IO.readlines(Gene).each{ |l|
				next if l =~ /^#/
				next if l =~ /^$/
				a = l.chomp.split("\t")
				raise("BED file does not have 6 columns. Aborting.") if a.size < 6
				sid, start, stop, name, score, strand = a
				outs1 << [sid, ".", Feature, start.to_i+1, stop, score, strand, "0", %|#{Attr} "#{name}"|]
				outs2 << [sid, start, stop, name, score, strand]
			}
			open("#{Adir}/#{Feature}.gtf", "w"){ |fout| 
				outs1.each{ |a| fout.puts a*"\t" }
			}
			open("#{Adir}/#{Feature}.bed", "w"){ |fout| 
				outs2.each{ |a| fout.puts a*"\t" }
			}
		end
	end
end
desc "01-5.make_sequence.gtf"
task "01-5.make_sequence.gtf", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

	outs1 = []
	outs2 = []
	IO.read(Fa).split(/^>/)[1..-1].each{ |ent|
		lab, *seq = ent.split("\n")
		name = lab.split(/\s+/)[0]
		len  = seq.join.gsub(/\s+/, "").size
		outs1 << [name, ".", "sequence", 1, len, 0, "+", "0", %|sequence_id "#{name}"|]
		outs2 << [name, 0, len, name, 0, "+"]
	}
	open("#{Adir}/sequence.gtf", "w"){ |fout| 
		outs1.each{ |a| fout.puts a*"\t" }
	}
	open("#{Adir}/sequence.bed", "w"){ |fout| 
		outs2.each{ |a| fout.puts a*"\t" }
	}
end
# }}} tasks 01


# {{{ tasks 02
desc "02-1.featureCounts"
task "02-1.featureCounts", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	mkdir_p Cdir

	Types.zip(Attrs){ |type, attr|
		anot = case Format 
					 when /^GFF3/ then "#{Adir}/#{type}.gff3"
					 else "#{Adir}/#{type}.gtf"
					 end
		bam  = "#{Mdir}/out.bam"
		fcnt = "#{Cdir}/#{type}.count"
		flog = "#{Cdir}/#{type}.count.log"
		cmd  = "featureCounts -T #{Threads} -t #{type} -g #{attr} -O -a #{anot} -o #{fcnt} #{bam} >#{flog} 2>&1"
		sh "#{cmd}"
	}
end
desc "02-2.bedcov"
task "02-2.bedcov", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

	Types.each{ |type|
		bed  = "#{Adir}/#{type}.bed"
		bam  = "#{Mdir}/out.bam"
		fcov = "#{Cdir}/#{type}.sumcov"
		flog = "#{Cdir}/#{type}.sumcov.log"
		cmd  = "samtools bedcov #{bed} #{bam} >#{fcov} 2>#{flog}"
		sh "#{cmd}"
	}
end
desc "02-3.normalize_and_make_output"
task "02-3.normalize_and_make_output", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

	script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	Types.each{ |type|
		fcnt = "#{Cdir}/#{type}.count"
		fcov = "#{Cdir}/#{type}.sumcov"
		frlt = "#{Cdir}/#{type}.result"
		cmd  = "ruby #{script} #{fcnt} #{fcov} #{frlt}"
		sh "#{cmd}"
	}
end
# }}} tasks 02


