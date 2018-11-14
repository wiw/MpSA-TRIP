<!DOCTYPE html>
<html>

<head>
    <meta charset="utf-8">
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

    <!-- Bootstrap CSS -->

    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/css/bootstrap.min.css" integrity="sha384-WskhaSGFgHYWDcbwN70/dfYBj47jz9qbsMId/iRN3ewGhXQFZCSftd1LZCfmhktB" crossorigin="anonymous">
    <title>TRIP $experiment report</title>
    <style>
        img {
            border: 2px;
            border-style: solid;
            border-color: #c3d9ff;
            padding: 5px;
        }
        
        h1,
        h2,
        h3,
        h4,
        h5,
        h6 {
            color: #000000;
            font-weight: bold;
        }
        
        p {
            text-indent: 2em;
            text-align: justify;
        }
        
        blockquote {
            background-color: #EEFFFF;
            border-left: 5px solid green;
            padding-left: 5px;
            margin: 0px;
            padding: 5px;
        }
        
        pre {
            background-color: #EEEEFF;
            display: block;
            border-left: 1px solid green;
            margin: 0px;
            padding: 5px;
        }
        
        code {
            background-color: #EEEEFF;
            margin: 15px;
            padding: 5px;
        }
        
        .long_string {
            word-break: break-all;
        }
        
        ol {
            counter-reset: item
        }
        
        ol > li {
            display: block
        }
        
        ol > li:before {
            content: counters(item, ".") " ";
            counter-increment: item
        }
        
        .suplement {
            color: gray;
        }
        
        .remark {
            font-size: 2em;
            font-weight: strong;
        }
    </style>
</head>

<body>
    <h1 style="text-align: center;">TRIP $experiment report</h1>
    <h2 style="text-align: center;" id="g0.1">$date / $lib</h2>
    <div class="container">
        <h3 id="g0.1.1">Used options:</h3>
        <blockquote class="long_string">
            $options
        </blockquote>
        <p></p>
        <h3 id="g0.1.2">Статистика количества ридов по всем индексам </h3>
        <p>${all_index_stat}</p>
        <p>
            <br>
        </p>
        <p></p>
        <h3 id="g0.1.3">I. Разделение исходных данных по индексам</h3>
        <p>Поскольку информация об illumina индексе есть только в прямых ридах, то выделить риды принадлежащие определенному индексу мы можем только установив соответствие между индексом и заголовками ридов. Для начала проводим анализ прямых ридов. Для этого каждую последовательность проверяем на соответствие следующим условиям:</p>
        <p>
        </p>
        <ul>
            <li>Полностью совпадающая индексная последовательность (0 ошибок в индексной последовательности) в начале рида</li>
            <li>Следующая за ней последовательность const1 с допустимыми 3 заменами</li>
        </ul>
        <div>Если находится рид подходящий под эти условия, то выполняются следующие действия:</div>
        <div>
            <ul>
                <li>Заголовок такого рида запоминается</li>
                <li>У рида отрезаем индекс и const1, получившаяся последовательность записывается в файл для дальнейшего анализа</li>
            </ul>
            <div>После того как все прямые риды будут обработаны, запускается анализ обратных ридов. Если заголовок обратного рида совпадает с заголовком, который находится в памяти, то такой обратный рид без изменений записывается в память.&nbsp;</div>
            <div>
                <br>
            </div>
            <div>
                ${current_index_stat}
                Общее количество ридов в файле: $TotalReadsCount</div>
            <div>
                <br>
            </div>
            <div>
                <br>
            </div>
        </div>
        <h3 id="g0.1.4">II. Выделение последовательностей баркода и генома из ридов</h3>
        <p>Для проведения анализа, необходимо ассоциировать баркод с геномной последовательностью. А также разделить данные по промоторным индексам и анализировать их независимо. Для этого было сделано следующее:</p>
        <p>
        </p>
        <p>Для каждого набора данных (прямого и обратного) выполняются следующие действия:</p>
        <ol>

            <li>Выбирается какое использовать регулярное выражение<b>*</b> для проверки ридов на основании анализируемого набора данных (прямого или обратного)</li>

            <li>Для каждого промоторного индекса из списка $pmi выполняются следующие действия:

                <ol>

                    <li>Создается пустой файл с меткой набора данных (прямой или обратный) и меткой обрабатываемого промоторного индекса для записи "процессированных" ридов</li>

                    <li>Для каждого рида из обрабатываемого набора данных:

                    </li>
                    <li>В память записывается результат работы регулярного выражения при обработке текущего рида</li>

                    <li>Если результат регулярного выражения не пустой, то:

                        <ol>

                            <li>Значащие части рида (если рид прямой, то значащие части: баркод, промоторный индекс, геном; если рид обратный, то значащие части: геном) в полученном результате проходят дополнительную проверку<b>**</b></li>

                            <li>Если дополнительная проверка прошла успешно тогда:

                                <ol>

                                    <li>Если рид прямой

                                        <ol>

                                            <li>Проверяется соответствует ли найденный промоторный индекс в риде промоторному индексу, который обрабатывается в настоящий момент и записывается в память под именем <b>matched_pmi</b>

                                                <ol>

                                                    <li>Если соответствует (или является мутантным вариантом), тогда в памяти создаются два пустых словаря следующей структуры: <b>А)</b> набор_данных -&gt; промоторный_индекс -&gt; заголовок_рида; <b>Б)</b> заголовок_рида -&gt; промоторный_индекс</li>

                                                    <li>Если не соответствует, то перейти к анализу следующего рида</li>
                                                </ol>
                                            </li>
                                        </ol>
                                    </li>

                                    <li>Если рид обратный, тогда:

                                        <ol>

                                            <li>Если заголовок рида существует в словаре <b>"Б"</b>, тогда из словаря <b>"Б"</b> выбирается промоторный индекс для этого заголовка и записывается в память под именем <b>matched_pmi</b>. А в словарь <b>"А"</b> для этого промоторного индекса записывается заголовок рида</li>

                                            <li>Если заголовка рида не существует в словаре <b>"Б"</b>, тогда перейти к анализу следующего рида</li>
                                        </ol>
                                    </li>

                                    <li>Если промоторный индекс <b>matched_pmi</b> из памяти соответствует индексу, который обрабатывается в настоящий момент, то:

                                        <ol>

                                            <li>Для каждой значащей части (кроме промоторного индекса, для прямых ридов: <i>barcode</i>, <i>genome</i>; для обратных: <i>genome</i>) выполняется следующее:

                                                <ol>

                                                    <li>В словарь <b>"А"</b> соответствующему набору данных, промоторному индексу и заголовку рида присваивается значащая часть</li>

                                                    <li>Отдельно, если значащая часть - это "<i>genome</i>", тогда последовательность генома в формате fastq записывается в файл</li>
                                                </ol>
                                            </li>
                                        </ol>
                                    </li>
                                </ol>
                            </li>

                            <li>Если дополнительная проверка не прошла, тогда перейти к анализу следующего рида</li>
                        </ol>
                    </li>

                </ol>
            </li>
        </ol>

        <div class="remark">*</div>
        <blockquote>
            <p>В зависимости от используемого набора данных, применяется различное регулярное выражение для поиска ридов. Условия по которым построены эти регулярные выражения написаны ниже.</p>
            &gt;
            <strong>Условия для поиска прямых ридов:</strong>
            <p>Значащие части рида: <i>barcode</i>, <i>pIndex</i>, <i>genome</i></p>
            <ul>

                <li>В начале рида есть неопределенное количество букв ("<i>barcode</i>"), после которых находится последовательность из четырех букв "<i>ATGC</i>" в любом сочетании и один любой индекс ("<i>pIndex</i>") из списка ${pmi_rev}. Количество возможных замен: ${pmiSubst}.</li>

                <li>Далее идет последовательность ${const_2}, допускается однонуклеотидная замена.</li>

                <li>Далее до конца рида идет геномная часть ("genome"), которая обязательно начинается с последовательности "GATC".</li>
            </ul>
            <strong>Условия для поиска обратных ридов:</strong>
            <p>Значащие части рида: genome</p>
            <ul>

                <li>В начале рида находится константная часть ${const_3}, количество возможных замен: ${const_2Error}.</li>

                <li>Далее до конца рида идет участок генома ("<i>genome</i>")</li>
            </ul>
        </blockquote>
        <div class="remark">**</div>
        <blockquote>
            <ul>

                <li>Проверка длины баркода, должна быть ${barcodeLength} +- ${barcodeError} нк.</li>

                <li>Проверка длины промоторного индекса, должна быть равна ${pmiLength}.</li>

                <li>В последовательности баркода и промоторного индекса не должно быть символов "N".</li>

                <li>Длина геномной части должна быть больше либо равна ${exp_fwd_genome_len}</li>
            </ul>
        </blockquote>

        <p></p>
        <hr>
        <div>Частотный анализ 4-х нуклеотидной последовательности:</div>
        <div>${spacer_4freq}</div>
        <hr>
        <p><span style="text-indent: 2em; font-size: 1rem; line-height: 1.5;">В вышеописанном алгоритме сначала обрабатываются прямые риды, а потом обратные. При этом не все прямые и обратные риды проходят установленные фильтры. Это приводит к тому, что возникает разница в количестве прямых и обратных ридов. Для парного выравнивания последовательностей на геном необходимо, чтобы количество прямых и обратных ридов совпадало и пары ридов были согласованы. Для этого заголовки прямых и обратных ридов снова загружаются в память и сравниваются между собой, оставляя только общие заголовки. Сформированный список общих заголовков используется устранения несогласованности между данными с прямыми и обратными ридами.</span>
            <br>
        </p>
        <table cellspacing="0" border="1px" bordercolor="#aaaacc" width="100%">
            <tbody>
                <tr>

                    <td>
                        <br>
                    </td>

                    <td>
                        прямые риды
                    </td>

                    <td>
                        обратные риды
                    </td>
                </tr>

                <tr>

                    <td>
                        до согласования
                    </td>

                    <td>
                        <br>
                    </td>

                    <td>
                        <br>
                    </td>
                </tr>
                <tr>

                    <td>
                        после согласования
                    </td>

                    <td>
                        <br>
                    </td>

                    <td>
                        <br>
                    </td>
                </tr>
            </tbody>
        </table>
        <p>
            Количество ридов потерянных на этом этапе:&nbsp;<span style="font-family: sans-serif;">${comm_difference}</span></p>
        <h3 id="g0.1.5">III. Выравнивание на геном</h3>
        <p>Для выравнивания на геном использовалась программа Bowtie2 со следующими параметрами:</p>
        <blockquote>
            <ul>
                <li><strong>Максимальное расстояние между ридами в паре</strong>&nbsp;- 5000 п.н.</li>
                <li><strong>Количество используемых ядер процессора</strong>&nbsp;- "максимально допустимое для данного процессора"</li>
                <li><strong>Поиск не более трех выравниваний</strong> для каждого рида <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#k-mode-search-for-one-or-more-alignments-report-each">-k mode (bowtie2 manual)</a></li>
                <li><strong>Не требуется полное совпадение</strong> рида при выравнивании <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#options">-local (Bowtie2 manual -&gt; Options -&gt; Alignment options)</a></li>
            </ul>
        </blockquote>
        <p><span style="text-indent: 2em; line-height: inherit; font-size: 1rem; text-align: left;">
<br>
</span></p>
        <p><span style="text-indent: 2em; line-height: inherit; font-size: 1rem; text-align: left;">Результаты выравнивания в свою очередь проходили многоуровневую фильтрацию:</span>
            <br>
        </p>
        <div>
            <ul>
                <li>Удаляются риды для которых не было найдено ни одного достоверного выравнивания</li>
                <li>Выбираются только те риды, который удовлетворяют следующим условиям:
                    <ul>
                        <li>для рида есть достоверное выравнивание</li>
                        <li>рид выравнивается как часть пары в парном выравнивании, а именно - если один рид из пары выравнивается на "+" цепь, то другой обязательно на "-" цепь</li>
                        <li>рид является частью пары и пара выравнивается согласованно, то есть находятся на одной хромосоме</li>
                    </ul>
                </li>

            </ul>
            <p><b>Statistic of alignment and filtering:</b></p>
            <table cellspacing="0" border="1px" bordercolor="#aaaacc" width="100%">
                <tbody>
                    <tr>
                        <td>&nbsp;Illumina index fst</td>
                        <td>alignment (результат работы Bowtie2)</td>
                        <td>filtering (сортировка и фильтрация после выравнивания)</td>
                    </tr>
                    <tr>
                        <td>promotor index fst</td>
                        <td>
                            <br>
                        </td>
                        <td>
                            <br>
                        </td>
                    </tr>
                    <tr>
                        <td>&gt;1 align&nbsp;
                            <br>
                        </td>
                        <td>
                            <br>
                        </td>
                        <td rowspan="3">
                            <br>
                        </td>
                    </tr>
                    <tr>
                        <td>1 align&nbsp;
                            <br>
                        </td>
                    </tr>
                    <tr>
                        <td>0 align&nbsp;
                            <br>
                        </td>
                    </tr>
                    <tr>
                        <td>promotor index snd</td>
                        <td>
                            <br>
                        </td>
                        <td>
                            <br>
                        </td>
                    </tr>
                    <tr>
                        <td>&gt;1 align</td>
                        <td>
                            <br>
                        </td>
                        <td rowspan="3">
                            <br>
                        </td>
                    </tr>
                    <tr>
                        <td>1 align</td>
                    </tr>
                    <tr>
                        <td>0 align</td>
                    </tr>
                    <tr>
                        <td><span style="font-family: sans-serif;">Illumina index snd</span>&nbsp;
                            <br>
                        </td>
                        <td>
                            <br>
                        </td>
                        <td>
                            <br>
                        </td>
                    </tr>
                    <tr>
                        <td>...etc</td>
                        <td>
                            <br>
                        </td>
                        <td>
                            <br>
                        </td>
                    </tr>
                </tbody>
            </table>
            <p><span style="text-indent: 2em; font-size: 1rem; line-height: 1.5;">
<br>
</span></p>
            <p><span style="text-indent: 2em; font-size: 1rem; line-height: 1.5;">Полученные данные конвертируются в текстовый файл ("bed"-формат, ${bed_file_path}), который затем загружается в память. Загруженный файл содержит в себе выровненные пары прямых и обратных ридов</span>
                <br>
            </p>
            <code class="long_string">
2L &nbsp; &nbsp; &nbsp;123862 &nbsp;123882 &nbsp;M02435:15:000000000-ARTHF:1:1101:21993:23024 &nbsp; &nbsp;255 &nbsp; &nbsp; + &nbsp; &nbsp; &nbsp; 99 &nbsp; &nbsp; &nbsp;20M &nbsp; &nbsp; = &nbsp; &nbsp; &nbsp; 123863 &nbsp;-62 &nbsp; &nbsp; GATCGTAATAGTCATAGATG &nbsp; &nbsp;GGGGGGGGGGGGGGGGGGGG &nbsp; &nbsp;AS:i:40 XN:i:0 &nbsp;XM:i:0 &nbsp;XO:i:0 &nbsp;XG:i:0 &nbsp;NM:i:0 &nbsp;MD:Z:20 YS:i:80 YT:Z:CP&nbsp;
</code></div>
        <div><code class="long_string">2L &nbsp; &nbsp; &nbsp;123862 &nbsp;123902 &nbsp;M02435:15:000000000-ARTHF:1:1101:21993:23024 &nbsp; &nbsp;255 &nbsp; &nbsp; - &nbsp; &nbsp; &nbsp; 147 &nbsp; &nbsp; 11S40M &nbsp;= &nbsp; &nbsp; &nbsp; 123863 &nbsp;62 &nbsp; &nbsp; &nbsp;CTCTTTGACTCGATCGTAATAGTCATAGATGTATACATGAGGAAAAGCATG &nbsp; &nbsp; FFFDEFGGGGGGGGGGGGGGGFGGGGGGFGGGGGGGGGGGGGGGGGGGGGG &nbsp; &nbsp; AS:i:80 XN:i:0 &nbsp;XM:i:0 &nbsp;XO:i:0 &nbsp;XG:i:0 &nbsp;NM:i:0 &nbsp;MD:Z:40 YS:i:40 YT:Z:CP</code>

            <p><span style="font-size: 1rem; line-height: 1.5; text-indent: 2em;">
<br>
</span></p>
            <p><span style="font-size: 1rem; line-height: 1.5; text-indent: 2em;">В рамках исследования нужная нам координата находится только в данных с обратными ридами, при этом геномная часть с как с прямых, так и с обратных ридов может выровнятся как на "+" так и на  "-" цепь. В выровненных данных необходимая координата располагается в начале (для  "+" цепи) или в конце (для "-" цепи) выравнивания длиной не менее 50 букв. В соответствии с этими условиями в память записывается словарь в котором каждому значению "Заголовок рида" (title) соответствует список из следующих пунктов: "Хромосома" (chr), "Координата" (coord), "Цепь" (strand).</span>
                <br>
            </p>
            <div>
                <h3 id="g0.1.6">IV. Сопоставление баркода и геномных данных</h3>
                <p>Для каждого промоторного индекса выполняются следующие действия:</p>

                <ul>
                    <li>Пересекаются данные по всем баркодам и выравниванию генома на основании общих заголовков ридов</li>
                    <li>
                        Производится анализ уникальных баркодов на наличие мутантных вариантов баркодов (****<i>Supl. 1</i>)
                    </li>
                    <li>На основании предыдущего пункта формируется словарь истинных баркодов, каждому из которых сопоставлен мутантный вариант этого баркода и геномная координата, которая ему соответствует</li>
                    <li>Далее на основании предыдущего пункта, для каждого истинного баркода производится поиск наиболее вероятной геномной координаты. Частота встречаемости, которой должна быть не ниже порога $mutationProbability от общей частоты встречаемости всех геномных координат в пределах этого истинного баркода</li>
                </ul>
                <p></p>
                <table cellspacing="0" border="1px" bordercolor="#aaaacc" width="100%">
                    <tbody>
                        <tr>

                            <td><span style="font-family: sans-serif;">Illumina index fst</span>
                                <br>
                            </td>

                            <td>reads count</td>

                            <td>effective reads count<b>***</b></td>

                            <td>unique barcode</td>

                            <td><span style="font-family: sans-serif;">unique barcode with genome alignment</span>
                                <br>
                            </td>
                        </tr>
                        <tr>

                            <td><span style="font-family: sans-serif;">Promotor index one</span>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>
                        </tr>
                        <tr>

                            <td><span style="font-family: sans-serif;">Promotor index two</span>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>
                        </tr>
                        <tr>

                            <td><span style="font-family: sans-serif;">Promotor index three</span>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>
                        </tr>
                        <tr>

                            <td>...etc.</td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>

                            <td>
                                <br>
                            </td>
                        </tr>
                    </tbody>
                </table>
                <blockquote>

                    <div class="remark">***</div>
                    Эффективное число ридов - число, которое является результатом сложения частот встречаемости всех истинных баркодов и частот встречаемости их мутантных вариантов. Отражает действительное число ридов использумых для анализа.</blockquote>
            </div>
            <h3 class="suplement" id="g0.1.7">****Supplementary material </h3>
            <h4 class="suplement" id="g0.1.7.1">Поиск мутантных вариантов баркода</h4>
            <div>
                <ul>
                    <li>Список всех баркодов преобразуется в сортированный по частоте встречаемости (от большего к меньшему) список с уникальными баркодами, при этом каждому повторяющемуся баркоду соответствует вычисленная частота встречаемости
                        <br>
                    </li>
                    <li>Далее оценивается длина большинства баркодов, если более чем 99% баркодов имеют одинаковую длину, то в дальнейшем анализируются только они. Иначе - анализируются все доступные.</li>
                    <li>Далее из списка баркодов убираются те у которых частота встречаемости меньше чем $readsValue</li>
                    <li>Затем если количество возможных ошибок в баркоде равно 0 - то весь список форматируется и записывается как есть. Если это не так, то работает алгоритм* поиска мутантных баркодов</li>
                    <li>*Для каждого баркода из данных сформированных на предыдущем шаге выполняется следующее:
                        <ul>
                            <li>Из последовательности баркода формируется список хешей (уникальные комбинации из символов "ATGC" характеризующие этот баркод или группу схожих баркодов)</li>
                            <li>Проверяется наличие хешей в списке уже обработанных баркодов. Если хоть один из хешей оказался в списке обработанных баркодов:

                                <ul>
                                    <li>Тогда для каждого такого хеша и соответствующего ему баркода производится поиск отличий (вычисление дистанции Левенштейна) между последовательностью баркода из хеша и анализируемой последовательностью баркода, если количество отличий не превышает $barcodeError - найденные баркоды записываются в список</li>
                                    <li>Если список возможных мутантных вариантов баркода непустой, тогда:
                                        <ul>
                                            <li>Если в этом списке более одного баркода, тогда в качестве мутантного выбираем баркод с максимальной частотой встречаемости.</li>
                                            <li>Если количество баркодов с максимальной частотой встречаемости более одного, то выбирается баркод, который выше в списке сортированных баркодов</li>
                                        </ul>
                                    </li>
                                    <li>Иначе - новые хеши записываются в список обработанных баркодов, а сам баркод записывается в словарь как уникальный
                                    </li>
                                </ul>
                            </li>
                            <li>Если таковых нет, то баркод считается - &nbsp;уникальным и немутированным. Хеши при этом добавляются в список обработанных баркодов в формате "хеш: баркод", а сам баркод записывается в словарь как уникальный.</li>
                        </ul>
                    </li>
                </ul>
            </div>

        </div>
    </div>
    <!-- Optional JavaScript -->
    <!-- jQuery first, then Popper.js, then Bootstrap JS -->
    <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js" integrity="sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49" crossorigin="anonymous"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/js/bootstrap.min.js" integrity="sha384-smHYKdLADwkXOn1EmN1qk/HfnUcbVRZyYmZ4qpPea6sjB/pTJ0euyQp0Mk8ck+5T" crossorigin="anonymous"></script>

</body>

</html>
