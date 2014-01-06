library(GEOmetadb)
con <- dbConnect(SQLite(), 'GEOmetadb.sqlite')

obQuery <- paste("SELECT DISTINCT gse.title,gse.gse",
                 "FROM",
                 "gse",
                 "WHERE",
                 "gse.type LIKE '%Expression profiling by array%' AND",
                 "gse.title LIKE '%obes%'", sep=" ")

obRes <- dbGetQuery(con, obQuery)


crQuery <- paste("SELECT DISTINCT gse.title,gse.gse",
                 "FROM",
                 "gse",
                 "WHERE",
                 "gse.type LIKE '%Expression profiling by array%' AND",
                 "gse.title LIKE '%crohn%'", sep=" ")

crRes <- dbGetQuery(con, crQuery)

metabQuery <- paste("SELECT DISTINCT gse.title,gse.gse",
                    "FROM",
                    "gse",
                    "WHERE",
                    "gse.type LIKE '%Expression profiling by array%' AND",
                    "gse.title LIKE '%metabol%'", sep=" ")

metabRes <- dbGetQuery(con, metabQuery)