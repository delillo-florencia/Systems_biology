expr <- expr |> 
  mutate(fc = case / control, 
         log2fc = log2(case / control))