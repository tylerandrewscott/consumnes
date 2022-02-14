https://www.webofscience.com/wos/woscc/summary/fb53728a-9312-44cc-a940-0b12ec99bf91-1ea7fc3b/relevance/1
install.packages("jsonlite")
install.packages("httr")
install.packages("base64enc")

library(devtools)
library(wosliterclient)

var.query.id <- 56 # integer | Retrieve records based on query identifier.
var.query.id <- 'fb53728a-9312-44cc-a940-0b12ec99bf91-1ea7fc3b'
var.count <- 10 # integer | Number of records to return, must be 0-100.
var.first.record <- 1 # integer | Specific record, if any within the result set to return. Cannot be less than 1 and greater than 100000.
var.sort.field <- 'sort.field_example' # character | Order by field(s). Field name and order by clause separated by '+', use A for ASC and D for DESC, ex: PY+D. Multiple values are separated by comma. If sortField was set on the original query, this parameter should be set as sorting is not a property of the query.

#Fetch record(s) by query identifier
api.instance <- SearchApi$new()
# Configure API key authorization: key
api.instance$apiClient$apiKeys['X-ApiKey'] <- 'clarivate:bcbc8640-3e7c-11ea-a5ca-bb71982fbd24';
result <- api.instance$QueryQueryIdGet(var.query.id, var.count, var.first.record, sort.field=var.sort.field)
dput(result)

test <- 'https://www.webofscience.com/wos/woscc/summary/fb53728a-9312-44cc-a940-0b12ec99bf91-1ea7fc3b/relevance/1'

library(rvest)
read_html(test)
wosliterclient::

https://wos-api.clarivate.com/api/wos
library(wosr)