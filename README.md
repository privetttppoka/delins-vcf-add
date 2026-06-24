# True delins / MNV detection

Набор скриптов для поиска истинных delins/MNV, которые в исходном VCF представлены как несколько соседних SNP. Такие случаи могут некорректно аннотироваться downstream-инструментами, например VEP, если оставить их отдельными SNP.

## Входные данные

Оба подхода используют:

- VCF или VCF.GZ с исходными вариантами;
- BAM или CRAM с выравниваниями;
- индекс для BAM/CRAM;
- индекс для VCF, если файл сжат и используется как indexed VCF;
- имя sample, если в VCF несколько образцов.

Скрипты анализируют только соседние biallelic SNP, прошедшие фильтр `PASS` и содержащие ALT-аллель в генотипе выбранного sample.

## Выходные файлы

Каждый скрипт формирует:

- TSV-файл с логом всех найденных кандидатов и решением по каждому из них;
- исправленный VCF/VCF.GZ, где подтвержденные группы соседних SNP заменены одной MNV-записью;
- опциональный JSON summary с общими счетчиками и временем выполнения.

Итоговый VCF остается нефазированным и подходит для дальнейшего downstream-анализа.

## Тестовые данные

В папке `indel_test` оставлен минимальный набор файлов для запуска примера:

- `indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.DeepVariant.vcf.gz` - входной VCF;
- `indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.DeepVariant.vcf.gz.csi` - индекс VCF;
- `indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.bam` - входной BAM;
- `indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.bam.bai` - индекс BAM.

Также добавлен пример результата после запуска `find_true_delins.py`:

- `example_true_delins.tsv` - лог найденных кандидатов;
- `example_true_delins.mnv.vcf.gz` - исправленный VCF с добавленными MNV;
- `example_true_delins.mnv.vcf.gz.tbi` - индекс итогового VCF.

Для сравнения добавлен пример результата после запуска `test_true_delins_whatshap.py`:

- `example_true_delins_whatshap.tsv` - лог кандидатов после анализа фазирования;
- `example_true_delins_whatshap.mnv.vcf.gz` - исправленный VCF по WhatsHap-подходу;
- `example_true_delins_whatshap.mnv.vcf.gz.tbi` - индекс итогового VCF.

## Скрипт `find_true_delins.py`

Основной вариант анализа через `pysam`.

Скрипт:

1. читает VCF и находит группы соседних SNP;
2. для каждого кандидата берет из BAM риды, перекрывающие все позиции SNP;
3. восстанавливает локальный гаплотип на каждом риде;
4. считает поддержку REF-, ALT- и смешанных гаплотипов;
5. добавляет подтвержденные delins/MNV в итоговый VCF.

Пример запуска:

```bash
python scripts/find_true_delins.py \
  --vcf indel_test/indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.DeepVariant.vcf.gz \
  --bam indel_test/indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.bam \
  --out-tsv indel_test/example_true_delins.tsv \
  --out-vcf indel_test/example_true_delins.mnv.vcf.gz
```

Аргументы:

- `--vcf` - входной VCF/VCF.GZ;
- `--bam` - входной BAM/CRAM с индексом;
- `--out-tsv` - таблица с кандидатами и статистикой;
- `--out-vcf` - исправленный VCF с добавленными MNV;
- `--sample` - имя sample, если нужно выбрать не первый sample в VCF;
- `--summary-json` - опциональный JSON с краткой сводкой.

## Скрипт `test_true_delins_whatshap.py`

Альтернативный вариант анализа через фазирование `WhatsHap`.

Скрипт:

1. временно приводит названия contig в VCF к стилю BAM, если они отличаются;
2. запускает `whatshap phase`;
3. анализирует фазированный VCF;
4. проверяет, находятся ли соседние SNP в одном phase set и лежат ли ALT-аллели на одном гаплотипе;
5. записывает итоговый нефазированный VCF с подтвержденными MNV.

Пример запуска:

```bash
python scripts/test_true_delins_whatshap.py \
  --vcf indel_test/indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.DeepVariant.vcf.gz \
  --bam indel_test/indel_test_chr11_118000000-119000000.MGI.cutadapt.bwa.MarkDuplicates.bam \
  --out-tsv indel_test/example_true_delins_whatshap.tsv \
  --out-vcf indel_test/example_true_delins_whatshap.mnv.vcf.gz
```

Аргументы:

- `--vcf` - входной VCF/VCF.GZ;
- `--bam` - BAM/CRAM для фазирования;
- `--out-tsv` - таблица с кандидатами и решением;
- `--out-vcf` - исправленный нефазированный VCF;
- `--sample` - имя sample, если нужно выбрать не первый sample в VCF;
- `--reference` - опциональный reference FASTA для WhatsHap;
- `--summary-json` - опциональный JSON с краткой сводкой.

Если `--reference` не указан, WhatsHap запускается с `--no-reference`.

## Общий модуль `true_delins_common.py`

Файл содержит общие функции:

- чтение и фильтрация SNP из VCF;
- работа с генотипами;
- построение кандидатов из соседних SNP;
- сопоставление названий contig между VCF и BAM;
- добавление INFO-полей и индексирование итогового VCF.

## Примечания

Скрипты автоматически обрабатывают разные стили названий хромосом, например `1` и `chr1`, если соответствие можно однозначно построить.

Подход через `pysam` обычно быстрее, потому что анализирует только области найденных кандидатов. Подход через `WhatsHap` медленнее, так как требует отдельного этапа фазирования, но может использоваться как альтернативная проверка.
