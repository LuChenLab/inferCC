package main


import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"log"
	"io"
	"os"
	"strconv"
	"strings"
	"sync"
)



func checkErr(err error) {
	if err != nil {
		panic(err)
	}
}


func readMeta(path string, groups []string) map[string][]string {
	log.Printf("Read meta %s", path)
	f, err := os.Open(path)
	defer f.Close()
	checkErr(err)

	res := make(map[string][]string)
	header := make(map[string]int)
	r := bufio.NewReader(f)

	for {
		line, err := r.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				checkErr(err)
			}
		}

		// clear csv format
		line = strings.Trim(line, "\n")
		line = strings.ReplaceAll(line, "\"", "")
		lines := strings.Split(line, ",")

		if len(header) == 0 {
			for i, j := range lines {
				header[j] = i
			}
		} else {
			// format cell name
			cellName := lines[header["cell_name"]]

			if cellName == "" {
				continue
			} else if cellName == "Macrophages" {
				cellName = "Monocytes"
			}

			lines[header["cell_name"]] = cellName

			// collect cell id by groups
			key := make([]string, 0, 0)
			for _, i := range groups {
				key = append(key, lines[header[i]])
			}

			temp, ok := res[strings.Join(key, "-")]
			if !ok {
				temp = []string{lines[0]}
			} else {
				temp = append(temp, lines[0])
			}
			res[strings.Join(key, "-")] = temp
		}
	}

	return res
}


func mean(data []float64) float64 {
	res := 0.0

	for _, i := range data {
		res += i
	}

	if len(data) > 0 {
		return res / float64(len(data))
	} else {
		return 0.0
	}
}


func calculateMean(wg *sync.WaitGroup, pipe chan string, output string, cellIdx map[string][]int) {
	defer wg.Done()
	log.Printf("Job %s start", output)
	f, err := os.OpenFile(output, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
	checkErr(err)
	defer f.Close()

	w := bufio.NewWriter(f)

	for {
		line, ok := <- pipe

		if !ok {
			break
		}

		line = strings.Trim(line, "\n")
		line = strings.ReplaceAll(line, "\"", "")
		lines := strings.Split(line, ",")
		
		newLine := []string{}
		for _, idx := range cellIdx {
			res := make([]float64, 0, 0)
			for _, j := range idx {
				temp, err := strconv.ParseFloat(lines[j], 32)
				checkErr(err)
				res = append(res, temp)
			}
			newLine = append(newLine, fmt.Sprintf("%v", mean(res)))
		}
		if len(newLine) != len(cellIdx) {
			log.Fatalf("%s %d not equal to header %d", lines[0], len(newLine), len(cellIdx))
		}
		w.WriteString(fmt.Sprintf("%s,%s\n", lines[0], strings.Join(newLine, ",")))
	}
	log.Printf("Job %s finished", output)
}



func readMain(infile, outfile string, nJobs int, group map[string][]string) {
	log.Printf("Read scale %s\n", infile)
	file, err := os.Open(infile)
	checkErr(err)
	gz, err := gzip.NewReader(file)
	checkErr(err)

	defer file.Close()
	defer gz.Close()

	r := bufio.NewReader(gz)

	wg := &sync.WaitGroup{}
	header := make(map[string][]int)

	pipes := make([]chan string, 0, 0)
	readLines := 0
	for {
		line, err := r.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				checkErr(err)
			}
		}

		if len(header) == 0 {
			line = strings.Trim(line, "\n")
			line = strings.ReplaceAll(line, "\"", "")
			lines := strings.Split(line, ",")
			
			// convert group - cell id to group - column index
			temp := make(map[string]int)
			for i, j := range lines {
				temp[j] = i
			}

			for i, j := range group {
				tempIdx := make([]int, 0, 0)
				for _, k := range j {
					tempIdx = append(tempIdx, temp[k])
				}
				header[i] = tempIdx
			}
			
			// init goroutines with different channel
			for i := 0; i < nJobs; i ++ {
				pipes = append(pipes, make(chan string, 1000))
				go calculateMean(wg, pipes[i], fmt.Sprintf("%s.%d", outfile, i), header)
				wg.Add(1)
			}
		} else {
			line = strings.TrimLeft(line, "gene\t")

			// read and process different row by different goroutines
			pipes[readLines % nJobs] <- line
			readLines ++

			if readLines % 100 == 0 {
				log.Printf("Read %d lines", readLines)
			}
		}
	}

	// close channel and wait task donw
	for _, i := range pipes {
		close(i)
	}
	wg.Wait()

	// collect all results to output file
	out, err := os.OpenFile(outfile, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
	checkErr(err)
	defer out.Close()
	w := bufio.NewWriter(out)

	temp_header := make([]string, 0, 0)
	for key := range header {
		temp_header = append(temp_header, key)
	}
	w.WriteString(fmt.Sprintf("gene,%s\n", strings.Join(temp_header, ",")))

	for i := 0; i < nJobs; i ++ {
		temp := fmt.Sprintf("%s.%d", outfile, i)

		r, err := os.Open(temp)
		checkErr(err)

		io.Copy(w, r)
		checkErr(os.Remove(temp))
	}
}


var (
	h bool
	meta string
	scale string
	output string
	group string
	nJobs int
)

func init() {
	flag.BoolVar(&h, "h", false, "this help")
	flag.StringVar(&meta, "m", "LungCancer10x/02_rds/meta_after_singleR.csv", "path to meta file")
	flag.StringVar(&scale, "s", "scale_data.csv.gz", "path to scaled data")
	flag.StringVar(&output, "o", "LungCancer10x/03_each_cells/mean_by_cell_disease_patient.tsv", "path output file")
	flag.StringVar(&group, "g", "cell_name,Disease,PatientID", "meta columns")
	flag.IntVar(&nJobs, "p", 3, "number of process")
}


func main() {
	flag.Parse()

    if h {
        flag.Usage()
    } else {
		cellIdx := readMeta(meta, strings.Split(group, ","))
	
		readMain(scale, output, nJobs, cellIdx)
	}
}