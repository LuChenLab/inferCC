package main

import (
	"crypto/tls"
	"fmt"
	"net/http"
	"os"
	"strings"

	"github.com/PuerkitoBio/goquery"
	"github.com/voxelbrain/goptions"
)

const targetUrl = "http://pms.cd120.com/download.html"

type Targets struct {
	Url      string
	Md5      string
	DataType string
	SampleID string
}

func checkErr(err error) {
	if err != nil {
		panic(err)
	}
}

func Get() []*Targets {

	tr := &http.Transport{
		TLSClientConfig: &tls.Config{InsecureSkipVerify: true},
	}

	client := &http.Client{Transport: tr}

	req, err := http.NewRequest("GET", targetUrl, nil)
	checkErr(err)

	resp, err := client.Do(req)
	checkErr(err)

	doc, err := goquery.NewDocumentFromReader(resp.Body)
	checkErr(err)

	dataTypes := make([]string, 10, 10)
	doc.Find("table").First().Find("th").Each(func(i int, selection *goquery.Selection) {
		dataTypes[i] = strings.Replace(strings.Trim(selection.Text(), " "), " ", "_", -1)
	})

	fmt.Print(fmt.Sprintf("%v\n", dataTypes))

	res := make([]*Targets, 0, 0)
	doc.Find("table").First().Find("tr").Each(func(i int, selection *goquery.Selection) {

	})

	return res
}

func main() {
	options := struct {
		Output string        `goptions:"-o, --output, description='Output path'"`
		Help   goptions.Help `goptions:"-h, --help, description='Show this help'"`
	}{}
	goptions.ParseAndFail(&options)

	if len(os.Args) <= 1 {
		goptions.PrintHelp()
		os.Exit(0)
	}

	Get()

}
